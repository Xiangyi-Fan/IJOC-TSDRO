yalmip clear;

timerVal = tic;
options = sdpsettings('solver', 'mosek');

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
kappa_hat = 1.9545e+04;
x_hat = [20.0000;
            52.9616;
            39.1447;
            100.0000;
            67.0773;
            52.9227;
            47.8937;
            100.0000];
    
% Parameters
T_hat = cell(ll,1);
T_hat{1} = [zeros(2*N,1);-1*kappa_hat];
T_hat{2} = zeros(aa,1);
G_hat = [eye(N),zeros(N),-1*x_hat;
         zeros(bb,aa)];
     
z_d = zeros(K,1);
pi_hat = cell(ll+1,K);
    
for k = 1:K
              
    %Solve by the solver
    % Decision Variables
    for l = 1:ll; Psi{l} = sdpvar(aa,aa);end
    M = sdpvar(aa,aa,'full');
    H = sdpvar(aa,aa,'full');
    Z = sdpvar(aa,aa,'full');
    
    % Constraints
    constraints2 = {};
    for l = 1:ll
        constraints2{end+1} = Psi{l} >= 0;
    end
        
    constraints2{end+1} = norm_inf(H) <= epsilon{k};
    constraints2{end+1} = norm_inf(Z) <= epsilon{k};
    cons1 = D*Omega{k}'+D*H'+F'*M;
    for l = 1:ll
        cons1 = cons1 - W{l}*Psi{l};
    end
    constraints2{end+1} = cons1 == 0;
        
    cons2 = eta{1}*Omega{k}'+eta{1}*Z';
    for l = 1:ll
        cons2 = cons2 - gamma{l,1}*Psi{l};
    end
    constraints2{end+1} = cons2 == 0;
    
    for l = 1:ll
        matrix1 = P{k}*Psi{l}*P{k}';
        constraints2{end+1} = matrix1(:) >= 0 ;
    end
    
    matrix2 = P{k}*M';
    constraints2{end+1} = matrix2(:) <= 0;
    constraints2{end+1} = M*ee <=0 ;
        
    % objective: sup
    cons3 = -1*trace(M'*G_hat);
    for l = 1:ll
        cons3 = cons3 + 0.5*trace(Psi{l}'*(ee*T_hat{l}'+T_hat{l}*ee'));
    end
    obj2 = cons3;
    obj2 = -1*obj2;
    % solving and post-processing
        
    out2 = optimize([constraints2{:}],obj2,options)

    z_d(k) = double(-1*obj2);
    for l = 1:ll
        pi_hat{l,k} = value(Psi{l});
    end
    pi_hat{ll+1,k} = value(M);
end
    
z_hat = kappa_hat + 1/delta*sum(z_d); 
   
optimal.obj = z_hat;
optimal.z = z_d;
optimal.t = toc(timerVal);