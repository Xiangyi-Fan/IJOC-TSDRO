function [x,obj] =  MS_PLD(data,epsilon0)

yalmip clear;

%% Parameters

N       = 8; % number of patients
Npoints = size(data,2); % number of data points
K = size(data,2); % number of partitions

l_l = 30; l_u = 60; % l_l <= l <= l_u
pi_l = 1 ; pi_u = 50; % pi_l <= pi <= pi_u 
c = 100;
T = (l_l+l_u)/2*N;

epsilon = cell(K,1);
for k = 1:K
    epsilon{k,1} = epsilon0/sqrt(Npoints);
end

% Define P_k
P = cell(K,1);
cons_points = data;
for k = 1:K
    P{k,1} = MS_PLD_partitions(cons_points,k);
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
for k = 1:K; Y{k} = sdpvar(a,b,'full');end
for k = 1:K; beta{k} = sdpvar(N,1);end
for k = 1:K; eta{k} = sdpvar(a,1);end
for k = 1:K; alpha{k} = sdpvar(p,N);end
for k = 1:K; gamma{k} = sdpvar(p,a);end
for k = 1:K; lambda{k} = sdpvar(1,1);end

%% Constraints

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = sum(x) <= T;

for k = 1:K

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
A = cell(K);
B = cell(K);
for k = 1:K
    A{k} = [zeros(N),eye(N),zeros(N,1)]'*[eye(N),zeros(N,1)]*Y{k};
    B{k} = c*[zeros(2*N,1);1]*[zeros(N,1);1]'*Y{k};
    value{k} = epsilon{k}*norm(A{k} + B{k} -lambda{k}*[zeros(2*N,1);1]*[zeros(2*N,1);1]','fro') + trace((A{k} + B{k} -lambda{k}*[zeros(2*N,1);1]*[zeros(2*N,1);1]')*Omega{k}) + lambda{k};
    values = values + value{k};
end
     
obj = values; % Infinite-norm + last component

%% solving and post-processing

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');

out = optimize([constraints{:}],obj,options);

double(obj)
double(x)

end