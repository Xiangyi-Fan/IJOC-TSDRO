function [optimal] = Dual_sp(cons_points,data,epsilon0,K,kappa_hat,x_hat)

yalmip clear;

timerVal = tic;
options = sdpsettings('verbose', 0, 'dualize',true, 'solver', 'mosek');

% Parameters

N       = 5; % number of stocks
Npoints = size(data,2); % number of data points

xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u
s_l = 0; s_u = 50; % s_l <= s <= s_u 
delta = 0.1; % risk attitude
B = 30; % total supply


g = [2.5;3;3.5;4;4.5]*2;

epsilon = cell(K,1);
for k = 1:K
    epsilon{k,1} = epsilon0/sqrt(Npoints);
end

% Define P_k
P = cell(K,1);
for k = 1:K
    P{k,1} = PLD_partitions(cons_points,k,K);
end

Omega = cell(K,1);
for k = 1:K
    I_k =  News_PLD_nums(P{k,1}, data);
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

vector = @(x) x(:); %transfer a matrix to a vector

% Index
aa = 2*N+1;
bb = N;
cc = 4*N;
ll = 2;
nn = 2;
mm = 1;


% Parameters
W = cell(ll,nn);
F = cell(nn,1);
gamma = cell(ll,mm);
W{1,1} = [zeros(bb,2*N),-1*g];
W{1,2} = [zeros(bb),-1*eye(bb),zeros(bb,1)];
W{2,1} = zeros(bb,aa);
W{2,2} = zeros(bb,aa);
gamma{1,1} = 1;
gamma{2,1} = 1;
F{1} = [eye(bb);zeros(bb);eye(bb);zeros(bb)];
F{2} = [zeros(bb);eye(bb);zeros(bb);eye(bb)];   
ee = [zeros(2*N,1);1];
T_hat = cell(ll,1);
T_hat{1} = [zeros(2*N,1);-1*kappa_hat];
T_hat{2} = zeros(aa,1);
G_hat = [-1*eye(bb),zeros(bb),x_hat;
         eye(bb),zeros(bb),-1*x_hat;
         zeros(bb,aa);
         zeros(bb,aa)];

     
z_d = zeros(K,1);
pi_hat = cell(ll+1,K);
runtime = zeros(K,1);
feasibility = zeros(K,1);

for k = 1:K
    
    Psi = cell(ll,1);
    
    % Decision Variables
    for l = 1:ll; Psi{l} = sdpvar(aa,aa);end
    M = sdpvar(cc,aa,'full');
    Z = sdpvar(aa,aa,'full');
    O = sdpvar(aa,aa);
    
    % Constraints
    constraints = {};
    
    constraints{end+1} = O >= 0;
    
    matrix2 = P{k}*O*ee;
    constraints{end+1} = matrix2(:) >= 0 ;    

    constraints{end+1} = Omega{k}+ Z - O == 0;

    constraints{end+1} =  ee'*O*ee == 1;
    
    for l = 1:ll
        constraints{end+1} = Psi{l} >= 0;
    end
        
    constraints{end+1} = norm(Z,'fro') <= epsilon{k};
    
    for n = 1:nn
        cons1 = F{n}'*M;
        for l = 1:ll
            cons1 = cons1 - W{l,n}*Psi{l};
        end
        constraints{end+1} = cons1 == 0;
    end

        
    cons2 = Omega{k}'+Z';
    for l = 1:ll
        cons2 = cons2 - gamma{l,1}*Psi{l};
    end
    constraints{end+1} = cons2 == 0;
    
    for l = 1:ll
        matrix1 = P{k}*Psi{l}*ee;
        constraints{end+1} = matrix1(:) >= 0 ;
    end
    
    matrix2 = P{k}*M';
    constraints{end+1} = matrix2(:) <= 0;
    constraints{end+1} = M*ee <=0 ;
        
    % objective: sup
    cons3 = -1*vector(M)'*vector(G_hat);
    for l = 1:ll
        cons3 = cons3 + 0.5*vector(Psi{l})'*vector(ee*T_hat{l}'+T_hat{l}*ee');
    end
    obj2 = cons3;
    obj2 = -1*obj2;
    
    % solving and post-processing
    out2 = optimize([constraints{:}],obj2,options);
    runtime(k) = out2.solvertime;
    feasibility(k) = out2.problem;

    z_d(k) = double(-1*obj2);
    pi_hat_k = cell(ll+1,1);
    for l = 1:ll
        pi_hat_k{l} = value(Psi{l});
    end
    pi_hat_k{ll+1} = value(M);
    pi_hat(:,k) = pi_hat_k; 
end
    
z_hat = kappa_hat;
for k = 1:K
    I_k =  News_PLD_nums(P{k,1}, data);
    z_hat = z_hat + 1/delta*z_d(k)*I_k/Npoints;
end
   
optimal.obj = z_hat;
optimal.z = z_d;
optimal.pi = pi_hat;
optimal.t = toc(timerVal);
optimal.runtime = max(runtime);
optimal.feasibility = feasibility;
end