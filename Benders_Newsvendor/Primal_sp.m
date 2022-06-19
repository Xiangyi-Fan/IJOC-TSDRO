function [optimal] = Primal_sp(data,epsilon0,K,kappa_hat,x_hat)

yalmip clear;

timerVal = tic;
options = sdpsettings('verbose', 0, 'dualize',false, 'solver', 'mosek');

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
cons_points = data;
for k = 1:K
    P{k,1} = PLD_partitions(cons_points,k,K);
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

vector = @(x) x(:); %transfer a matrix to a vector

% Index
aa = 2*N+1;
bb = N;
cc = 4*N;
dd = size(P{k,1},1);
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

% Parameters
ee = [zeros(2*N,1);1];
    
% Parameters
T_hat = cell(ll,1);
T_hat{1} = [zeros(2*N,1);-1*kappa_hat];
T_hat{2} = zeros(aa,1);
G_hat = [-1*eye(bb),zeros(bb),x_hat;
         eye(bb),zeros(bb),-1*x_hat;
         zeros(bb,aa);
         zeros(bb,aa)];
     
z_p = zeros(K,1);
    
for k = 1:K
              
    % Decision Variables
    for l = 1:ll; theta{l} = sdpvar(1);end
    for l = 1:ll; sigma{l} = sdpvar(dd,1);end
    alpha = sdpvar(dd,cc,'full');
    beta = sdpvar(cc,1);
    Y1 = sdpvar(bb,aa,'full');
    Y2 = sdpvar(bb,aa,'full');
    Q = sdpvar(aa,aa,'full');
    
    % Constraints
    constraints1 = {};
    for l = 1:ll
        constraints1{end+1} = theta{l} >= 0;
        constraints1{end+1} = sigma{l}(:) >= 0;
    end
    
    constraints1{end+1} = alpha(:) >= 0;
    constraints1{end+1} = beta >= 0;
    
    for l = 1:ll
        constraints1{end+1} = 0.5*((W{l,1}'*Y1+W{l,2}'*Y2+gamma{l,1}*Q-ee*T_hat{l}')+(W{l,1}'*Y1+W{l,2}'*Y2+gamma{l,1}*Q-ee*T_hat{l}')')-theta{l}*ee*ee'-0.5*(P{k,1}'*sigma{l}*ee'+ee*sigma{l}'*P{k,1}) >= 0;
    end
    
     
    constraints1{end+1} = F{1}*Y1+F{2}*Y2-G_hat-alpha'*P{k,1}-beta*ee' == 0;
        
    % objective: sup
    obj1 = trace(Q*Omega{k,1})+epsilon{k,1}*norm_1(Q);
    
    
    % solving and post-processing
        
    out1 = optimize([constraints1{:}],obj1,options);

    z_p(k) = double(obj1);

end
    
z_hat = kappa_hat + 1/delta*sum(z_p); 
   
optimal.obj = z_hat;
optimal.z = z_p;
optimal.t = toc(timerVal);

end