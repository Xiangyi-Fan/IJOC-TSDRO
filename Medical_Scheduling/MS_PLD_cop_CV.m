function [x,obj] =  MS_PLD_cop_CV(cons_points,data,epsilon0_cop,K)

yalmip clear;

%% Parameters

N       = 8; % number of patients
Npoints = size(data,2); % number of data points

l_l = 20; l_u = 100; % l_l <= l <= l_u
c = 200;
T = (l_l+l_u)/2*N;
delta = 0.1; % risk attitude


% Define P_k
P = cell(K,1);
for k = 1:K
    P{k,1} = MS_PLD_partitions(cons_points,k,K);
end

% K fold CV
epsilon = cell(K,1);
for k = 1:K
    I_k = MS_PLD_nums(P{k,1}, data);
    if I_k > 0
        epsilon{k,1} = epsilon0_cop/(sqrt(I_k)*K);
    else
        epsilon{k,1} = epsilon0_cop;
    end
end


Omega = cell(K,1);
num = cell(K,1);
for k = 1:K
    Omega{k,1} = zeros(size(data,1),size(data,1));
    num{k} = 0;
    for i = 1:Npoints
        if P{k,1}*data(:,i) >= 0
            Omega{k,1} = Omega{k,1} + data(:,i)*data(:,i)';
            num{k} = num{k} + 1;
        end
    end
    if num{k} > 0
        Omega{k,1} = Omega{k,1}/num{k};
    end
end


% Matrix D
D = [];
for i = 1:N
    vec = zeros(N+1,1);
    vec(i) = -1;
    vec(i+1) = 1;
    D = [D;vec'];
end

vector = @(x) x(:); %transfer a matrix to a vector

%% Decision Variables

p = size(P{1},1);
a = N+1;
b = 2*N+1;
x = sdpvar(N,1);
kappa = sdpvar(1,1);
for k = 1:K; Q{k} = sdpvar(2*N+1,2*N+1,'full');end
for k = 1:K; theta{k} = sdpvar(1,1,'full');end
for k = 1:K; Y{k} = sdpvar(a,b,'full');end
for k = 1:K; phi{k} = sdpvar(1,1,'full');end
for k = 1:K; V1{k} = sdpvar(b,b,'full');end
for k = 1:K; V2{k} = sdpvar(b,b,'full');end
for k = 1:K; V3{k} = sdpvar(b,b,'full');end
for k = 1:K; W1{k} = sdpvar(b,b);end
for k = 1:K; W2{k} = sdpvar(b,b);end
for k = 1:K; W3{k} = sdpvar(b,b);end
for k = 1:K; Sigma1{k} = sdpvar(p,p);end
for k = 1:K; Sigma2{k} = sdpvar(p,p);end
for k = 1:K; Sigma3{k} = sdpvar(p,p);end
for k = 1:K; beta{k} = sdpvar(N,1);end
for k = 1:K; eta{k} = sdpvar(a,1);end
for k = 1:K; alpha{k} = sdpvar(p,N,'full');end
for k = 1:K; gamma{k} = sdpvar(p,a,'full');end
for k = 1:K; lambda{k} = sdpvar(1,1);end
for k = 1:K; B{k} = sdpvar(2*N+1,2*N+1);end
for k = 1:K; pi{k} = sdpvar(1,1);end


%% Constraints

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = sum(x) <= T;

A1 = cell(K);
B1 = cell(K);

for k = 1:K
    %% 
    
    A1{k} = [zeros(N),eye(N),zeros(N,1)]'*[eye(N),zeros(N,1)]*Y{k};
    B1{k} = c*[zeros(2*N,1);1]*[zeros(N,1);1]'*Y{k};
    constraints{end+1} = theta{k}+kappa >= 0; 
    constraints{end+1} = 0.5*((Q{k}-A1{k}-B1{k})+(Q{k}-A1{k}-B1{k})')-theta{k}*[zeros(2*N,1);1]*[zeros(2*N,1);1]' == V2{k};

    %% 

    constraints{end+1} = phi{k} >= 0;
    constraints{end+1} = 0.5*(Q{k}+Q{k}')-phi{k}*[zeros(2*N,1);1]*[zeros(2*N,1);1]' == V3{k};

    
    %% 
    
    constraints{end+1} = B{k}+pi{k}*[zeros(2*N,1);1]*[zeros(2*N,1);1]' == V1{k};
    
    %% SDP Approximation

    constraints{end+1} = V1{k} == W1{k}+P{k}'*Sigma1{k}*P{k};
    constraints{end+1} = Sigma1{k}(:) >= 0;
    constraints{end+1} = W1{k} >= 0;
    constraints{end+1} = V2{k} == W2{k}+P{k}'*Sigma2{k}*P{k};
    constraints{end+1} = Sigma2{k}(:) >= 0;
    constraints{end+1} = W2{k} >= 0;
    constraints{end+1} = V3{k} == W3{k}+P{k}'*Sigma3{k}*P{k};
    constraints{end+1} = Sigma3{k}(:) >= 0;
    constraints{end+1} = W3{k} >= 0;


    %% 

    constraints{end+1} = alpha{k} >= 0;
    constraints{end+1} = beta{k}+x >= 0;
    constraints{end+1} = D*Y{k} - [eye(N),zeros(N,N+1)] - alpha{k}'*P{k} - beta{k}*[zeros(2*N,1);1]' == 0;

    %%
    constraints{end+1} = gamma{k} >= 0;
    constraints{end+1} = eta{k} >= 0;
    constraints{end+1} = Y{k} - gamma{k}'*P{k} - eta{k}*[zeros(2*N,1);1]' == 0;

end

%% objective: min

value = cell(K,1);
values = 0;
for k = 1:K
    value{k} = pi{k}+epsilon{k}*norm(Q{k}'+B{k},'fro') + vector(Q{k}')'*vector(Omega{k})+ vector(B{k})'*vector(Omega{k});
    values = values + num{k}*value{k};
end
     
obj = kappa+1/delta*1/Npoints*values; % Infinite-norm + last component

%% solving and post-processing

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');

out = optimize([constraints{:}],obj,options);
runtime = out.solvertime;

end