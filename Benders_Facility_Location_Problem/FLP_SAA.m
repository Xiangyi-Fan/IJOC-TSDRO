function [x_value,obj_value,runtime] =  FLP_SAA(data,f,C)

yalmip clear;

options = sdpsettings('verbose', 0, 'dualize',false, 'solver', 'mosek');

% Parameters

I = 5;
J = 5;
Npoints = size(data,2); % number of data points

u = 160*ones(I,1);
g = 2000*ones(J,1);


% Matrix D
DD = [];
for i = 1:I
    matx = [];
    for j = 1:J
        ee0 = zeros(J+1,1);
        ee0(j) = C(i,j);
        matx = [matx;ee0'];
    end
    DD = [DD;matx];
end
matx = [];
for j = 1:J
    ee0 = zeros(J+1,1);
    ee0(j) = g(j);
    matx = [matx;ee0'];
end
DD = [DD;matx];

% Index
ll = I+2*J+I*J;
N_1 = Npoints;
N_2 = I*J+J;

% decision variables
x = binvar(I,1);
Y = sdpvar(N_2,N_1,'full');
Y_ = reshape(Y,N_1*N_2,1);

% Parameters
W_hat = cell(ll,1);
T_hat = cell(ll,1);
for l = 1:J
    T_hat{l} = ee(J+1,J+1);
    W_hat{l} = [];
    for i = 1:I+1
        W_hat{l} = [W_hat{l};[zeros(J),ee(J,l)]];
    end
end

for l = J+1:J+I
    l1 = l - J;
    T_hat{l} = -1*ee(J+1,J+1)*u(l1)*x(l1);
    W_hat{l} = [];
    for i = 1:I+1
        if i ~= l1
            W_hat{l} = [W_hat{l};zeros(J,J+1)];
        else
            W_hat{l} = [W_hat{l};[-1*eye(J),zeros(J,1)]];
        end
    end
end
for l = J+I+1:2*J+I+I*J
    l1 = l - (J+I);
    T_hat{l} = zeros(J+1,1);
    W_hat{l} = zeros(I*J+J,J+1);
    W_hat{l}(l1,J+1) = 1;
end 

  
xi = data(1:J+1,:);


constraints = {};

r0 = DD';
T0 = [];
for l = 1:ll
    T0 = [T0;T_hat{l}'];
end

R_ = kron(speye(N_1),r0);
T_ = kron(speye(N_1),T0);
xi_ = reshape(xi,N_1*(J+1),1);
Q0 = [];
for i=1:N_1
    Q = sparse(zeros(1,N_1));
    Q(1,i) = 1;
    W0 = zeros(ll,N_2);
    for l = 1:ll
        W0(l,:) = (W_hat{l}*xi(:,i))';
    end
    Q = kron(Q,W0);
    Q0 = [Q0;Q];
end
W_ = Q0;

obj0 = f'*x + 1/N_1*xi_'*R_*Y_;
constraints{end+1} = T_*xi_ <= W_*Y_;


obj = obj0;

out = optimize([constraints{:}],obj,options);
runtime = out.solvertime;

x_value = double(x);
obj_value = double(obj);

end