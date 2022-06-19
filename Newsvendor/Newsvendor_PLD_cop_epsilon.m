function [x,obj,runtime] =  Newsvendor_PLD_cop_epsilon(cons_points,data,epsilon0_cop,K)

yalmip clear;

%% Parameters

N       = 5; % number of stocks
Npoints = size(data,2); % number of data points

xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u
s_l = 0; s_u = 50; % s_l <= s <= s_u 
delta = 0.1; % risk attitude
B1 = 30; % total supply
epsilon1 = 0.05; % mean & convariance ambiguity set 

R = 0.1; % sqrt(N*xi_u*xi_u + N*s_u*s_u +1)
delta2 = 0.1;
epsilon2 = R^2/sqrt(Npoints)*(2+sqrt(2*log(1/delta2))); % Frobenius norm ambiguity set

b0 = max(xi_u,s_u);
b0 = max(1,b0);
a0 = min(xi_l,s_l);
a0 = min(1,a0);
delta3 = 0.1;
epsilon3 = 0.01*(b0-a0)*sqrt((log((2*N)^2/2+2*N)-log(delta3))/(2*Npoints)); % Infinite norm ambiguity set  

g = [2.5;3;3.5;4;4.5]*2;

% Define P_k
P = cell(K,1);
for k = 1:K
    P{k,1} = PLD_partitions(cons_points,k,K);
end

epsilon = cell(K,1);
for k = 1:K
    I_k = News_PLD_nums(P{k,1}, data);
    R_k = News_PLD_radius(P{k,1}, cons_points, k);
    if I_k > 0
        epsilon{k,1} = epsilon0_cop{k}*R_k^2/sqrt(I_k);
    else
        epsilon{k,1} = inf;
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

vector = @(x) x(:); %transfer a matrix to a vector

%% Decision Variables

l = size(P{1},1);
x = sdpvar(N,1);
kappa = sdpvar(1,1);
for k = 1:K; Q{k} = sdpvar(2*N+1,2*N+1,'full');end
for k = 1:K; theta{k} = sdpvar(1,1,'full');end
for k = 1:K; Y1{k} = sdpvar(N,2*N+1,'full');end
for k = 1:K; Y2{k} = sdpvar(N,2*N+1,'full');end
for k = 1:K; beta{k} = sdpvar(1,1,'full');end
for k = 1:K; V1{k} = sdpvar(2*N+1,2*N+1,'full');end
for k = 1:K; V2{k} = sdpvar(2*N+1,2*N+1,'full');end
for k = 1:K; V3{k} = sdpvar(2*N+1,2*N+1,'full');end
for k = 1:K; W1{k} = sdpvar(2*N+1,2*N+1);end
for k = 1:K; W2{k} = sdpvar(2*N+1,2*N+1);end
for k = 1:K; W3{k} = sdpvar(2*N+1,2*N+1);end
for k = 1:K; Sigma1{k} = sdpvar(l,l);end
for k = 1:K; Sigma2{k} = sdpvar(l,l);end
for k = 1:K; Sigma3{k} = sdpvar(l,l);end
for k = 1:K; n{k} = sdpvar(N,1);end
for k = 1:K; b{k} = sdpvar(N,1);end
for k = 1:K; d{k} = sdpvar(N,1);end
for k = 1:K; h{k} = sdpvar(N,1);end
for k = 1:K; M{k} = sdpvar(l,N);end
for k = 1:K; A{k} = sdpvar(l,N);end
for k = 1:K; C{k} = sdpvar(l,N);end
for k = 1:K; F{k} = sdpvar(l,N);end
for k = 1:K; pi{k} = sdpvar(1,1);end
for k = 1:K; B{k} = sdpvar(2*N+1,2*N+1);end

%% Constraints

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = sum(x) <= B1;

for k = 1:K
    %% 

    constraints{end+1} = theta{k}+kappa >= 0; 
    constraints{end+1} = 0.5*((Q{k}-[zeros(2*N,1);1]*g'*Y1{k}-[zeros(N);eye(N);zeros(1,N)]*Y2{k})+(Q{k}-[zeros(2*N,1);1]*g'*Y1{k}-[zeros(N);eye(N);zeros(1,N)]*Y2{k})')-theta{k}*[zeros(2*N,1);1]*[zeros(2*N,1);1]' == V2{k};

    %% 

    constraints{end+1} = beta{k} >= 0;
    constraints{end+1} = 0.5*(Q{k}+Q{k}')-beta{k}*[zeros(2*N,1);1]*[zeros(2*N,1);1]' == V3{k};
    
    %% 
    
    constraints{end+1} = B{k} + pi{k}*[zeros(2*N,1);1]*[zeros(2*N,1);1]' == V1{k};
   
    
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

    constraints{end+1} = M{k} >= 0;
    constraints{end+1} = n{k}-x >= 0;
    constraints{end+1} = Y1{k}+[eye(N),zeros(N,N+1)]-M{k}'*P{k}-n{k}*[zeros(2*N,1);1]' == 0;

    %%
    constraints{end+1} = A{k} >= 0;
    constraints{end+1} = b{k}+x >= 0;
    constraints{end+1} = Y2{k}-[eye(N),zeros(N,N+1)]-A{k}'*P{k}-b{k}*[zeros(2*N,1);1]' == 0;

    %%

    constraints{end+1} = C{k} >= 0;
    constraints{end+1} = d{k} >= 0;
    constraints{end+1} = Y1{k}-C{k}'*P{k}-d{k}*[zeros(2*N,1);1]' == 0;

    %%

    constraints{end+1} = F{k} >= 0;
    constraints{end+1} = h{k} >= 0;
    constraints{end+1} = Y2{k}-F{k}'*P{k}-h{k}*[zeros(2*N,1);1]' == 0;
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