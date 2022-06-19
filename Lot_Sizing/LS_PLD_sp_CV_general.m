function [x,obj] =  LS_PLD_sp_CV_general(cons_points,data,epsilon0)

yalmip clear;

%% Parameters

N       = 5; % number of locations
Npoints = size(data,2); % number of data points
K = size(data,2); % number of partitions

v_l = 40; v_u = 50; % v_l <= v <= v_u
u_l = 20 ; u_u = 40; % u_l <= u <= u_u 
C = [40;50;60;70;80];
cap = 80;


% Define P_k
P = cell(K,1);
for k = 1:K
    P{k,1} = LS_PLD_partitions(cons_points,k);
end

% K fold CV
epsilon = cell(K,1);
for k = 1:K
    I_k = LS_PLD_nums(P{k,1}, data);
    if I_k > 0
        epsilon{k,1} = epsilon0/(sqrt(I_k)*K);
    else
        epsilon{k,1} = 9999999999999;
    end
end

epsilon_p = 0;


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


% Matrix A and B
A1 = [];
B1 = [];
for i = 1:N
    vec=zeros(N,1);
    positions1=i;
    vec(positions1)=1;
    A0 = repmat(vec,N,1);
    A1 = [A1;A0'];
    B0 = zeros(N*N,1);
    positions2 = [(i-1)*N+1:i*N];
    B0(positions2)=-1;
    B1 = [B1;B0'];
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
for k = 1:K; V1{k} = sdpvar(b,b,'full');end
for k = 1:K; W1{k} = sdpvar(b,b);end
for k = 1:K; theta1{k} = sdpvar(l,1);end
for k = 1:K; B{k} = sdpvar(b,b);end
for k = 1:K; pi{k} = sdpvar(1,1);end
eta0 = sdpvar(1,1);
omega = sdpvar(1,1);
for k = 1:K; r{k} = sdpvar(1,1);end
for k = 1:K; s{k} = sdpvar(1,1);end


%% Constraints

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = x <= cap;
constraints{end+1} = omega >= 0;

for k = 1:K

    %% 

    constraints{end+1} = alpha{k} >= 0;
    constraints{end+1} = beta{k}+x >= 0;
    constraints{end+1} = A1*Y{k} + B1*Y{k} - [zeros(N,N*N),eye(N),zeros(N,1)] - alpha{k}'*P{k} - beta{k}*[zeros(a+N,1);1]' == 0;

    %%
    constraints{end+1} = gamma{k} >= 0;
    constraints{end+1} = eta{k} >= 0;
    constraints{end+1} = Y{k} - gamma{k}'*P{k} - eta{k}*[zeros(a+N,1);1]' == 0;

    %% 
    
    constraints{end+1} = B{k} + pi{k}*[zeros(N+N*N,1);1]*[zeros(N+N*N,1);1]' == V1{k};
    constraints{end+1} = V1{k} == W1{k}+0.5*(P{k}'*theta1{k}*[zeros(N+N*N,1);1]'+[zeros(N+N*N,1);1]*theta1{k}'*P{k});
    constraints{end+1} = theta1{k} >= 0;
    constraints{end+1} = W1{k} >= 0;       

    %%
    constraints{end+1} = s{k} + eta0 <= omega;
    constraints{end+1} = sqrt(4*r{k}^2+(s{k}+eta0)^2) <= (2*omega-s{k}-eta0);
    value = cell(K,1);
    value{k} = pi{k}+epsilon{k}*norm([eye(a);zeros(N+1,a)]*Y{k}+B{k},'fro') + vector(([eye(a);zeros(N+1,a)]*Y{k})')'*vector(Omega{k})+ vector(B{k})'*vector(Omega{k});
    constraints{end+1} = value{k} <= s{k};
    
end

%% objective: min
     
obj = C'*x + epsilon_p*omega - eta0 + 2*omega*1; % Infinite-norm + last component
for k = 1:K
    obj = obj - 2/Npoints*num{k}*r{k};
end

%% solving and post-processing

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');

out = optimize([constraints{:}],obj,options);

end