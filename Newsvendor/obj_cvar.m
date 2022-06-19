function [obj] =  obj_cvar(x,data,g,Spoints)

N = size(x,1);
delta = 0.1; % risk attitude

beta = sdpvar(1,1);

constraints = {};


N_1 = Spoints;
xi0 = data(1:2*N,:);
r0 = [1;zeros(N,1);zeros(N,1)];
h0 = [0;-beta;x;-x;zeros(N,1);zeros(N,1)];
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

obj = beta + 1/delta*obj0;

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');
out = optimize([constraints{:}],obj,options);

end

