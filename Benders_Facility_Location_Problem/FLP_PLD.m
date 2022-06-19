function [x,obj,runtime] =  LS_PLD(data,epsilon0)

yalmip clear;

%% Parameters

N       = 5; % number of locations
Npoints = size(data,2); % number of data points
K = size(data,2); % number of partitions

v_l = 40; v_u = 50; % v_l <= v <= v_u
u_l = 20 ; u_u = 40; % u_l <= u <= u_u 
C = [50;60;70;80;90]*0.8;
cap = 80;

epsilon = cell(K,1);
for k = 1:K
    epsilon{k,1} = epsilon0/sqrt(Npoints);
end

% Define P_k
P = cell(K,1);
cons_points = data;
for k = 1:K
    P{k,1} = LS_PLD_partitions(cons_points,k);
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

% Matrix A and B
A = [];
B = [];
for i = 1:N
    vec=zeros(N,1);
    positions1=i;
    vec(positions1)=1;
    A0 = repmat(vec,N,1);
    A = [A;A0'];
    B0 = zeros(N*N,1);
    positions2 = [(i-1)*N+1:i*N];
    B0(positions2)=-1;
    B = [B;B0'];
end

vector = @(x) x(:); %transfer a matrix to a vector

%% Decision Variables

l = size(P{1},1);
a = N*N;
b = (N*N+N+1);
x = sdpvar(N,1);
for k = 1:K; Y{k} = sdpvar(a,b,'full');end
for k = 1:K; beta{k} = sdpvar(N,1);end
for k = 1:K; eta{k} = sdpvar(a,1);end
for k = 1:K; alpha{k} = sdpvar(l,N);end
for k = 1:K; gamma{k} = sdpvar(l,a);end

%% Constraints

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = x <= cap;

for k = 1:K

    %% 

    constraints{end+1} = alpha{k} >= 0;
    constraints{end+1} = beta{k}+x >= 0;
    constraints{end+1} = A*Y{k} + B*Y{k} - [zeros(N,N*N),eye(N),zeros(N,1)] - alpha{k}'*P{k} - beta{k}*[zeros(a+N,1);1]' == 0;

    %%
    constraints{end+1} = gamma{k} >= 0;
    constraints{end+1} = eta{k} >= 0;
    constraints{end+1} = Y{k} - gamma{k}'*P{k} - eta{k}*[zeros(a+N,1);1]' == 0;
    
end

%% objective: min

value = cell(K,1);
values = 0;
for k = 1:K
    value{k} = epsilon{k}*norm_1([eye(a);zeros(N+1,a)]*Y{k}) + vector(([eye(a);zeros(N+1,a)]*Y{k})')'*vector(Omega{k});
    values = values + value{k};
end
     
obj = C'*x + values; % Infinite-norm + last component

%% solving and post-processing

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');

out = optimize([constraints{:}],obj,options);
runtime = out.solvertime;

end