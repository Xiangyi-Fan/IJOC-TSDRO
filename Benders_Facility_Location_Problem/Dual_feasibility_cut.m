function [cut] = Dual_feasibility_cut(cons_points,data,epsilon0,k,x_hat,C)

yalmip clear;

options = sdpsettings('verbose', 0, 'dualize',true, 'solver', 'mosek');

% Parameters

I = 5;
J = 5;
Npoints = size(data,2); % number of data points

delta = 1; % risk attitude


u = 160*ones(I,1);
g = 2000*ones(J,1);


epsilon = epsilon0/sqrt(Npoints);

% Define P_k
P = FLP_PLD_partitions(cons_points,k);


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
     
H = cell(ll,1);
    
% Decision Variables
for l = 1:ll; H{l} = sdpvar(J+1,J+1);end
G = sdpvar(J+1,J+1,'full');
O = sdpvar(J+1,J+1);

% Constraints
constraints2 = {};

constraints2{end+1} = O >= 0;
    
matrix2 = P*O*ee(J+1,J+1);
constraints2{end+1} = matrix2(:) >= 0 ;    

constraints2{end+1} = G - O == 0;

constraints2{end+1} =  ee(J+1,J+1)'*O*ee(J+1,J+1) == 0;
    
for l = 1:ll
    constraints2{end+1} = H{l} >= 0;
end

constraints2{end+1} = norm(G,'fro') <= epsilon;
    
for l = 1:ll
    constraints2{end+1} =  ee(J+1,J+1)'*H{l}*ee(J+1,J+1) >= 0;
end
    
cons2 = DD*G';
for l = 1:ll
    cons2 = cons2 - W_hat{l}*H{l};
end
constraints2{end+1} = cons2 == 0;
    
    
for l = 1:ll
    matrix1 = P*H{l}*ee(J+1,J+1);
    constraints2{end+1} = matrix1(:) >= 0 ;
end


for l = 1:ll
    constraints2{end+1} = H{l}(:) >= -1*ones((J+1)^2,1) ;
    constraints2{end+1} = H{l}(:) <= 1*ones((J+1)^2,1);
end
constraints2{end+1} = G(:) >= -1*ones((J+1)^2,1) ;
constraints2{end+1} = G(:) <= 1*ones((J+1)^2,1) ;
constraints2{end+1} = O(:) >= -1*ones((J+1)^2,1) ;
constraints2{end+1} = O(:) <= 1*ones((J+1)^2,1) ;

        
% objective: sup
cons3 = 0;
for l = 1:ll
    cons3 = cons3 + 0.5*vector(H{l})'*vector(ee(J+1,J+1)*T_hat{l}'+T_hat{l}*ee(J+1,J+1)');
end
obj2 = cons3;
obj2 = -1*obj2;
    
% solving and post-processing
out2 = optimize([constraints2{:}],obj2,options);

z_d = double(-1*obj2)
pi_hat = value(H); 

   
cut.pi = pi_hat;
cut.runtime = out2.solvertime;
end