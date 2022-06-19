function [x,obj] =  LS_PLD_SAA(data)

yalmip clear;

%% Parameters

N       = 5; % number of stocks
Npoints = size(data,2); % number of data points
K = Npoints;
a = N*N;

C = [50;60;70;80;90]*0.8;
cap = 80;

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
  
v = data(1:N*N,:);
u = data(N*N+1:N*N+N,:);

x = sdpvar(N,1);

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = x <= cap;

% Define P_k
P = cell(K,1);
cons_points = data;
for k = 1:K
    P{k,1} = LS_PLD_partitions(cons_points,k);
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
    r0 = [eye(a);zeros(N+1,a)];
    h0 = [-1*x];
    T0 = [zeros(N,a),eye(N),zeros(N,1)];
    W0 = [A+B];
    M0 = size(T0,1);
    d0 = a+N;
    d1 = a+N+1;
    N_2 = size(r0,2);
    R_ = kron(speye(N_1),r0);
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
for k = 1:K
    value{k} = trace([eye(a);zeros(N+1,a)]*Y{k}*Omega{k});
    values = values + value{k};
end

obj0 = C'*x + values;
obj = obj0;

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');
out = optimize([constraints{:}],obj,options);

double(x)
double(obj)
end