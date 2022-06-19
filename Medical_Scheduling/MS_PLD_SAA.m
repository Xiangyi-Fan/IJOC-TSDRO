function [x,obj] =  MS_PLD_SAA(data)

yalmip clear;

%% Parameters

N       = 8; % number of patients
Npoints = size(data,2); % number of data points
K = Npoints;

l_l = 30; l_u = 60; % l_l <= l <= l_u
c = 100;
T = (l_l+l_u)/2*N;

% Matrix D
D = [];
for i = 1:N
    vec = zeros(N+1,1);
    vec(i) = -1;
    vec(i+1) = 1;
    D = [D;vec'];
end
  

x = sdpvar(N,1);

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = sum(x) <= T;

% Define P_k
P = cell(K,1);
cons_points = data;
for k = 1:K
    P{k,1} = MS_PLD_partitions(cons_points,k);
end

Y = cell(K,1);
Y_ = cell(K,1);
Xi = cell(K,1);
Omega = cell(K,1);
for k = 1:K
    Omega{k,1} = zeros(size(data,1),size(data,1));
    Xi{k,1} = [];
    for i = 1:Npoints
        if P{k,1}*data(:,i) >= 0
            Omega{k,1} = Omega{k,1} + data(:,i)*data(:,i)';
            Xi{k,1} = [Xi{k,1},data(:,i)];
        end
    end
    Omega{k,1} = Omega{k,1}/Npoints;
end



for k = 1:K
    xi = Xi{k,1};
    xi0 = xi;
    N_1 = size(xi,2);
    r0 = [zeros(N),eye(N),zeros(N,1)]'*[eye(N),zeros(N,1)];
    h0 = [-1*x];
    T0 = [eye(N),zeros(N,N+1)];
    W0 = D;
    M0 = size(T0,1);
    d1 = N+N+1;
    N_2 = size(r0,2);
    T_ = kron(speye(N_1),T0);
    xi_ = reshape(xi0,N_1*d1,1);
    h_ = repmat(h0,N_1,1);
    W_ = kron(speye(N_1),W0);
    Y{k} = sdpvar(N_2,d1,'full');
    Y_{k} = reshape(Y{k}*xi,N_1*N_2,1);

    constraints{end+1} = T_*xi_+h_ <= W_*Y_{k};
    constraints{end+1} = Y_{k} >= 0;
end



value = cell(K,1);
values = 0;
A = cell(K);
B = cell(K);
for k = 1:K
    A{k} = [zeros(N),eye(N),zeros(N,1)]'*[eye(N),zeros(N,1)]*Y{k};
    B{k} = c*[zeros(2*N,1);1]*[zeros(N,1);1]'*Y{k};
    value{k} = trace((A{k} + B{k})*Omega{k});
    values = values + value{k};
end
obj = values;



options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');
out = optimize([constraints{:}],obj,options)

double(x)
double(obj)
end