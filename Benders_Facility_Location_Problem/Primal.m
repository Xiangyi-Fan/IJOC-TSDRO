function [optimal] = Primal(data,epsilon0,K,x_hat,C)

yalmip clear;

timerVal = tic;
options = sdpsettings('verbose', 0, 'dualize',false, 'solver', 'mosek');

% Parameters

I = 5;
J = 5;
Npoints = size(data,2); % number of data points

delta = 1; % risk attitude


u = 160*ones(I,1);
g = 2000*ones(J,1);


epsilon = cell(K,1);
for l = 1:K
    epsilon{l,1} = epsilon0/sqrt(Npoints);
end

% Define P_k
P = cell(K,1);
cons_points = data;
for l = 1:K
    P{l,1} = FLP_PLD_partitions(cons_points,l);
end

Omega = cell(K,1);
for l = 1:K
    Omega{l,1} = zeros(size(data,1),size(data,1));
    for i = 1:Npoints
        if P{l,1}*data(:,i) >= 0
            Omega{l,1} = Omega{l,1} + data(:,i)*data(:,i)';
        end
    end
end

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


vector = @(x) x(:); %transfer a matrix to a vector

% Index
ll = I+2*J+I*J;


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
    T_hat{l} = -1*ee(J+1,J+1)*u(l1)*x_hat(l1);
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


p0 = size(P{1},1);
z_p = zeros(K,1);
pi_hat = cell(ll,K);
runtime = zeros(K,1);
feasibility = zeros(K,1);

for k = 1:K
    
    
    % Decision Variables
    for l = 1:ll; pi{l} = sdpvar(1,1);end
    Y = sdpvar(I*J+J,J+1);
    B = sdpvar(J+1,J+1);
    gamma = sdpvar(1);
    V1 = sdpvar(J+1,J+1);
    for l = 1:ll; V2{l} = sdpvar(J+1,J+1);end
    W1 = sdpvar(J+1,J+1);
    for l = 1:ll; W2{l} = sdpvar(J+1,J+1);end
    theta1 = sdpvar(p0,1);
    for l = 1:ll; theta2{l} = sdpvar(p0,1);end
  
    % Constraints
    constraints2 = {};

    for l = 1:ll
        constraints2{end+1} = pi{l} >= 0;
    end
    
    constraints2{end+1} = B + gamma*ee(J+1,J+1)*ee(J+1,J+1)' == V1;
    
    for l = 1:ll
        constraints2{end+1} = 0.5*(W_hat{l}'*Y + Y'*W_hat{l}-ee(J+1,J+1)*T_hat{l}'-T_hat{l}*ee(J+1,J+1)') - pi{l}*ee(J+1,J+1)*ee(J+1,J+1)' == V2{l};
    end
    
    constraints2{end+1} = V1 == W1+0.5*(P{k}'*theta1*ee(J+1,J+1)'+ee(J+1,J+1)*theta1'*P{k});
    constraints2{end+1} = theta1 >= 0;
    constraints2{end+1} = W1 >= 0; 
    
    for l = 1:ll
        constraints2{end+1} = V2{l} == W2{l}+0.5*(P{k}'*theta2{l}*ee(J+1,J+1)'+ee(J+1,J+1)*theta2{l}'*P{k});
        constraints2{end+1} = theta2{l} >= 0;
        constraints2{end+1} = W2{l} >= 0;
    end

        
    % objective: sup
    
    values = gamma+epsilon{k}*norm(Y'*DD+B,'fro') + vector(Y'*DD)'*vector(Omega{k})+ vector(B)'*vector(Omega{k});
     
    obj2 = 1/Npoints*values; % Infinite-norm + last component
    
    % solving and post-processing
    out2 = optimize([constraints2{:}],obj2,options);
    runtime(k) = out2.solvertime;
    feasibility(k) = out2.problem;
    z_p(k) = double(obj2);
end    

z_hat = sum(z_p); 
   
optimal.obj = z_hat;
optimal.z = z_p;
optimal.t = toc(timerVal);
optimal.runtime = max(runtime);
optimal.feasibility = feasibility;
end