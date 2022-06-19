function [optimal] = Primal(data,epsilon0,K,kappa_hat,x_hat)

yalmip clear;

timerVal = tic;
options = sdpsettings('verbose', 0, 'dualize',false,'solver', 'mosek');

%% Parameters

N       = 8; % number of patients
Npoints = size(data,2); % number of data points

l_l = 20; l_u = 100; % l_l <= l <= l_u
c = 200;
TT = (l_l+l_u)/2*N;
delta = 0.1; % risk attitude


epsilon = cell(K,1);
for k = 1:K
    epsilon{k,1} = epsilon0/sqrt(Npoints);
end

% Define P_k
P = cell(K,1);
cons_points = data;
for k = 1:K
    P{k,1} = MS_PLD_partitions(cons_points,k,K);
end

Omega = cell(K,1);
for k = 1:K
    Omega{k,1} = zeros(size(data,1),size(data,1));
    for i = 1:Npoints
        if P{k,1}*data(:,i) >= 0
            Omega{k,1} = Omega{k,1} + data(:,i)*data(:,i)';
        end
    end
    Omega{k,1} = Omega{k,1}/Npoints;
end

% Matrix D
DD = [];
for i = 1:N
    vec = zeros(N+1,1);
    vec(i) = -1;
    vec(i+1) = 1;
    DD = [DD;vec'];
end

A = [zeros(N),eye(N),zeros(N,1)]'*[eye(N),zeros(N,1)];
B = c*[zeros(2*N,1);1]*[zeros(N,1);1]';

%% Index
aa = 2*N+1;
bb = N+1;
cc = size(P{1,1},1);
ll = 2;
nn = 1;
mm = 1;

%% step 0

% Parameters
W = cell(ll,1);
gamma = cell(ll,nn);
eta = cell(mm,1);
D = zeros(bb,aa);
eta{1} = 1;
W{1} = -1*(A+B)';
W{2} = zeros(bb,aa);
gamma{1,1} = 1;
gamma{2,1} = 1;
F = [DD; eye(bb)];
    
% Parameters
ee = [zeros(2*N,1);1];
    
% Parameters
T_hat = cell(ll,1);
T_hat{1} = [zeros(2*N,1);-1*kappa_hat];
T_hat{2} = zeros(aa,1);
G_hat = [eye(N),zeros(N),-1*x_hat;
         zeros(bb,aa)];
     
z_p = zeros(K,1);
pi_hat = cell(ll+1,K);
    
for k = 1:K
              
    %Solve by the solver
    % Decision Variables
    for l = 1:ll; theta{l} = sdpvar(1);end
    for l = 1:ll; Sigma{l} = sdpvar(cc,cc);end
    alpha = sdpvar(cc,aa,'full');
    beta = sdpvar(aa,1);
    Y = sdpvar(bb,aa,'full');
    Q = sdpvar(aa,aa,'full');
    
    % Constraints
    constraints1 = {};
    for l = 1:ll
        constraints1{end+1} = theta{l} >= 0;
        constraints1{end+1} = Sigma{l}(:) >= 0;
    end
    
    constraints1{end+1} = alpha(:) >= 0;
    constraints1{end+1} = beta >= 0;
    
    for l = 1:ll
        constraints1{end+1} = 0.5*((W{l}'*Y+gamma{l,1}*Q-ee*T_hat{l}')+(W{l}'*Y+gamma{l,1}*Q-ee*T_hat{l}')')-theta{l}*ee*ee'-P{k,1}'*Sigma{l}*P{k,1} >= 0;
    end
    
     
    constraints1{end+1} = F*Y-G_hat-alpha'*P{k,1}-beta*ee' == 0;
        
    % objective: sup
    obj1 = trace(D'*Y*Omega{k,1})+epsilon{k,1}*norm_1(D'*Y)+eta{1,1}*(trace(Q*Omega{k,1})+epsilon{k,1}*norm_1(Q));
    
    
    % solving and post-processing
        
    out1 = optimize([constraints1{:}],obj1,options);

    z_p(k) = double(obj1);

end
    
z_hat = kappa_hat + 1/delta*sum(z_p); 
   
optimal.obj = z_hat;
optimal.z = z_p;
optimal.t = toc(timerVal);

end