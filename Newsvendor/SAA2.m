function [x,obj,runtime] =  SAA2(data)

yalmip clear;

%% Parameters

N       = 5; % number of stocks
Npoints = size(data,2); % number of data points

B = 30; % total supply
g = [2.5;3;3.5;4;4.5]*2;
  
xi = data(1:N,:);
s = data(N+1:2*N,:);
delta = 0.1; % risk attitude

kappa = sdpvar(1,1);
x = sdpvar(N,1);

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = sum(x) <= B; 


N_1 = Npoints;
xi0 = [xi;s];
r0 = [1;zeros(N,1);zeros(N,1)];
h0 = [0;-kappa;x;-x;zeros(N,1);zeros(N,1)];
T0 = [zeros(1,2*N);zeros(1,2*N);[-eye(N);zeros(N)]';[eye(N);zeros(N)]';zeros(N,2*N);zeros(N,2*N)];
W0 = [[1,zeros(1,N),zeros(1,N)];[1,-g',zeros(1,N)];[zeros(N,1),eye(N),zeros(N)];[zeros(N,1),zeros(N),eye(N)];[zeros(N,1),eye(N),zeros(N)];[zeros(N,1),zeros(N),eye(N)]];
M0 = size(T0,1);
d0 = 2*N;
N_2 = size(r0,1);
Q0 = [];
obj0 = 0;
T_ = kron(speye(N_1),T0);
xi_ = reshape(xi0,N_1*d0,1);
h_ = repmat(h0,N_1,1);
W_ = kron(speye(N_1),W0);
r_ = repmat(r0,N_1,1);
Y = sdpvar(N_2,N_1,'full');
Y_ = reshape(Y,N_1*N_2,1);

for i=1:N_1
    Q = sparse(zeros(1,N_1));
    Q(1,i) = 1;
    Q = kron(Q,[zeros(1,N_2);[zeros(N+1,1);-[zeros(N);eye(N)]'*xi0(:,i)]';zeros(2*N,N_2);zeros(2*N,N_2)]);
    Q0 = [Q0;Q];
end
W_ = W_ + Q0;
obj0 = 1/N_1*r_'*Y_;
constraints{end+1} = T_*xi_+h_ <= W_*Y_;

obj = kappa + 1/delta*obj0;

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');
out = optimize([constraints{:}],obj,options);
runtime = out.solvertime;

end 