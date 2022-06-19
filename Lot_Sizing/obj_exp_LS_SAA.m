function [obj,feasibility] =  obj_exp_LS_SAA(x,data)

yalmip clear;

%% Parameters

N       = 5; % number of stocks
Spoints = size(data,2); % number of data points
a = N*N;

C = [40;50;60;70;80];

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

constraints = {};

N_1 = Spoints;
xi0 = [v;u];
r0 = [eye(a);zeros(N,a)];
h0 = [-1*x];
T0 = [zeros(N,a),eye(N)];
W0 = [A+B];
M0 = size(T0,1);
d0 = a+N;
N_2 = size(r0,2);
R_ = kron(speye(N_1),r0);
T_ = kron(speye(N_1),T0);
xi_ = reshape(xi0,N_1*d0,1);
h_ = repmat(h0,N_1,1);
W_ = kron(speye(N_1),W0);
Y = sdpvar(N_2,N_1,'full');
Y_ = reshape(Y,N_1*N_2,1);

obj0 = C'*x + 1/N_1*xi_'*R_*Y_;
constraints{end+1} = T_*xi_+h_ <= W_*Y_;
constraints{end+1} = Y_ >= 0;

obj = obj0;

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');
out = optimize([constraints{:}],obj,options);
feasibility = out.problem;


end