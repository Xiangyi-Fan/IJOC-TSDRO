function [optimal] = Dual(cons_points,data,epsilon0,K,x_hat,C)

yalmip clear;

timerVal = tic;
options = sdpsettings('verbose', 0, 'dualize',true, 'solver', 'mosek');

% Parameters

I = 5;
J = 5;
Npoints = size(data,2); % number of data points

delta = 1; % risk attitude

u = 160*ones(I,1);
g = 2000*ones(J,1);

% K fold CV
epsilon = cell(K,1);
for k = 1:K
    epsilon{k,1} = epsilon0/sqrt(Npoints);
end

% Define P_k
P = cell(K,1);
for k = 1:K
    P{k,1} = FLP_PLD_partitions(cons_points,k);
end

Omega = cell(K,1);
for k = 1:K
    I_k =  FLP_PLD_nums(P{k,1}, data);
    Omega{k,1} = zeros(size(data,1),size(data,1));
    for i = 1:Npoints
        if P{k,1}*data(:,i) >= 0
            Omega{k,1} = Omega{k,1} + data(:,i)*data(:,i)';
        end
    end
    if I_k > 0
        Omega{k,1} = Omega{k,1}/I_k;
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
lambda_hat = cell(ll,1);
kappa_hat = cell(ll,1);
W_hat = cell(ll,1);
T_hat = cell(ll,1);
for l = 1:ll-1
    lambda_hat{l} = 0;
    kappa_hat{l} = 0;
end
lambda_hat{ll} = 1;

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


z_d = zeros(K,1);
pi_hat = cell(ll,K);
runtime = zeros(K,1);
feasibility = zeros(K,1);

for k = 1:K
    
    H = cell(ll,1);
    
    % Decision Variables
    for l = 1:ll; H{l} = sdpvar(J+1,J+1);end
    G = sdpvar(J+1,J+1,'full');
    O = sdpvar(J+1,J+1);
    
    % Constraints
    constraints2 = {};
    
    constraints2{end+1} = O >= 0;
    
    matrix2 = P{k}*O*ee(J+1,J+1);
    constraints2{end+1} = matrix2(:) >= 0 ;    

    constraints2{end+1} = Omega{k}+ G - O == 0;

    constraints2{end+1} =  ee(J+1,J+1)'*O*ee(J+1,J+1) == 1;
    
    for l = 1:ll
        constraints2{end+1} = H{l} >= 0;
    end
        
    constraints2{end+1} = norm(G,'fro') <= epsilon{k};
    
    for l = 1:ll
        constraints2{end+1} =  ee(J+1,J+1)'*H{l}*ee(J+1,J+1) >= 0;
    end
        
    cons2 = DD*Omega{k}'+ DD*G';
    for l = 1:ll
        cons2 = cons2 - W_hat{l}*H{l};
    end
    constraints2{end+1} = cons2 == 0;
    
    for l = 1:ll
        matrix1 = P{k}*H{l}*ee(J+1,J+1);
        constraints2{end+1} = matrix1(:) >= 0 ;
    end
  
        
    % objective: sup
    cons3 = 0;
    for l = 1:ll
        cons3 = cons3 + 0.5*vector(H{l})'*vector(ee(J+1,J+1)*T_hat{l}'+T_hat{l}*ee(J+1,J+1)');
    end
    obj2 = cons3;
    obj2 = -1*obj2;
    
    % solving and post-processing
    out2 = optimize([constraints2{:}],obj2,options);
    runtime(k) = out2.solvertime;
    feasibility(k) = out2.problem;

    z_d(k) = double(-1*obj2);
    pi_hat_k = cell(ll,1);
    for l = 1:ll
        pi_hat_k{l} = value(H{l});
    end
    pi_hat(:,k) = pi_hat_k; 
end
    
z_hat = 0;
for k = 1:K
    I_k =  FLP_PLD_nums(P{k,1}, data);
    z_hat = z_hat + I_k/Npoints*z_d(k);   
end
   
optimal.obj = z_hat;
optimal.z = z_d;
optimal.pi = pi_hat;
optimal.t = toc(timerVal);
optimal.runtime = max(runtime);
optimal.feasibility = feasibility;
end