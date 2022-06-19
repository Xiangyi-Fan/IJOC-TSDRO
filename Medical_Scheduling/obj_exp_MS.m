function [obj] =  obj_exp_MS(x,data)

yalmip clear;

%% Parameters

N       = 8; % number of stocks
Spoints = size(data,2); % number of data points

c = 100;


% Matrix D
D = [];
for i = 1:N
    vec = zeros(N+1,1);
    vec(i) = -1;
    vec(i+1) = 1;
    D = [D;vec'];
end
  

constraints = {};

N_1 = Spoints;
xi0 = data;
r0 = [zeros(N),eye(N),zeros(N,1)]'*[eye(N),zeros(N,1)];
r1 = c*[zeros(N,1);1];
h0 = [-1*x];
T0 = [eye(N),zeros(N,N+1)];
W0 = D;
M0 = size(T0,1);
d0 = N+N+1;
N_2 = size(r0,2);
R0_ = kron(speye(N_1),r0);
R1_ = repmat(r1,N_1,1);
T_ = kron(speye(N_1),T0);
xi_ = reshape(xi0,N_1*d0,1);
h_ = repmat(h0,N_1,1);
W_ = kron(speye(N_1),W0);
Y = sdpvar(N_2,N_1,'full');
Y_ = reshape(Y,N_1*N_2,1);

obj0 = 1/N_1*(xi_'*R0_*Y_ + R1_'*Y_);
constraints{end+1} = T_*xi_+h_ <= W_*Y_;
constraints{end+1} = Y_ >= 0;

obj = obj0;

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');
out = optimize([constraints{:}],obj,options);

double(obj)

end