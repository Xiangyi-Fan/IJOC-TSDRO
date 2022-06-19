yalmip clear;
timerVal = tic;

%% Parameters

N       = 8; % number of patients
Npoints = size(data,2); % number of data points

l_l = 20; l_u = 100; % l_l <= l <= l_u
c = 200;
T = (l_l+l_u)/2*N;
delta = 0.1; % risk attitude


epsilon = cell(K,1);
for k = 1:K
    epsilon{k,1} = epsilon0_cop/sqrt(Npoints);
end

% Define P_k
P = cell(K,1);
cons_points = data;
for k = 1:K
    P{k,1} = MS_PLD_partitions(cons_points,k,K);
end

Omega = cell(K,1);
for k = 1:K
    Omega{k,1} = zeros(size(data,1),size(data,1));
    for i = 1:Npoints
        if P{k,1}*data(:,i) >= 0
            Omega{k,1} = Omega{k,1} + data(:,i)*data(:,i)';
        end
    end
    Omega{k,1} = Omega{k,1}/Npoints;
end

% Matrix D
D = [];
for i = 1:N
    vec = zeros(N+1,1);
    vec(i) = -1;
    vec(i+1) = 1;
    D = [D;vec'];
end

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
for k = 1:K; V2{k} = sdpvar(b,b,'full');end
for k = 1:K; V3{k} = sdpvar(b,b,'full');end
for k = 1:K; W2{k} = sdpvar(b,b);end
for k = 1:K; W3{k} = sdpvar(b,b);end
for k = 1:K; Sigma2{k} = sdpvar(p,p);end
for k = 1:K; Sigma3{k} = sdpvar(p,p);end
for k = 1:K; beta{k} = sdpvar(N,1);end
for k = 1:K; eta{k} = sdpvar(a,1);end
for k = 1:K; alpha{k} = sdpvar(p,N,'full');end
for k = 1:K; gamma{k} = sdpvar(p,a,'full');end
for k = 1:K; lambda{k} = sdpvar(1,1);end

%% Constraints

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = sum(x) <= T;

A = cell(K);
B = cell(K);

for k = 1:K
    %% 
    
    A{k} = [zeros(N),eye(N),zeros(N,1)]'*[eye(N),zeros(N,1)]*Y{k};
    B{k} = c*[zeros(2*N,1);1]*[zeros(N,1);1]'*Y{k};
    constraints{end+1} = theta{k}+kappa >= 0; 
    constraints{end+1} = 0.5*((Q{k}-A{k}-B{k})+(Q{k}-A{k}-B{k})')-theta{k}*[zeros(2*N,1);1]*[zeros(2*N,1);1]' == V2{k};

    %% 

    constraints{end+1} = phi{k} >= 0;
    constraints{end+1} = 0.5*(Q{k}+Q{k}')-phi{k}*[zeros(2*N,1);1]*[zeros(2*N,1);1]' == V3{k};

    %% SDP Approximation

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
%     value{k} = epsilon{k}*norm(Q{k}-lambda{k}*[zeros(2*N,1);1]*[zeros(2*N,1);1]',1) + trace((Q{k}-lambda{k}*[zeros(2*N,1);1]*[zeros(2*N,1);1]')*Omega{k}) + lambda{k};
    value{k} = epsilon{k}*norm(Q{k},1) + trace(Q{k}*Omega{k});
    values = values + value{k};
end
     
obj = kappa+1/delta*values; % Infinite-norm + last component

%% solving and post-processing

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');

out = optimize([constraints{:}],obj,options)

double(kappa)
double(obj)
double(x)
toc(timerVal)