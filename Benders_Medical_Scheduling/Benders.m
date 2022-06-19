function [optimal] = Benders(cons_points,data,epsilon0,K)

yalmip clear;

timerVal = tic;
options = sdpsettings('verbose', 0, 'solver', 'mosek');


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
for k = 1:K
    P{k,1} = MS_PLD_partitions(cons_points,k,K);
end

Omega = cell(K,1);
for k = 1:K
    I_k =  MS_PLD_nums(P{k,1}, data);
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
z_u = Inf;
x_opt = zeros(N,1);
U = zeros(aa,aa,ll+1,0,K);
% V = zeros(aa,aa,ll+1,0,K);
obj_list = [];
error = 0.05;
count = zeros(1,K);

% Parameters
T = cell(ll,1);
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



%% begin
runtime = 0;
runtime_sub = 0;
logical = true;
iter = 0;
while logical
    
    %% step 1
    
    % Decision Variables
    x = sdpvar(N,1);
    r = sdpvar(1,1,K);
    kappa = sdpvar(1,1);
    
    % Parameters
    T{1} = [zeros(2*N,1);-1*kappa];
    T{2} = zeros(aa,1);
    G = [eye(N),zeros(N),-1*x;
         zeros(bb,aa)];
    ee = [zeros(2*N,1);1];
     
    
     
    % Constraints
    constraints1 = {};
    
    constraints1{end+1} = kappa >= -9999;
    for k = 1:K
        constraints1{end+1} = r(1,1,k) >= -9999;
    end
    
    constraints1{end+1} = x >= 0;
    constraints1{end+1} = sum(x) <= TT;
    
    for k = 1:K
        for j = 1:size(U,4)
            if ~isnan(U(:,:,1,j,k))
                cons = -1*trace(G'*U(:,:,3,j,k));
                for l = 1:ll
                    cons = cons + 0.5*trace(U(:,:,l,j,k)'*(ee*T{l}'+T{l}*ee'));
                end
                constraints1{end+1} = cons <= r(:,:,k);
            end
        end  
    end

    % objective: min
    obj = kappa;
    for k = 1:K
        I_k =  MS_PLD_nums(P{k,1}, data);
        obj = obj + 1/delta*I_k/Npoints*r(:,:,k);   
    end
    
    % solving and post-processing
    out = optimize([constraints1{:}],obj,options);
    runtime = runtime + out.solvertime;
    z_l = double(obj);
    x_hat = value(x);
    kappa_hat = value(kappa);
    r_hat = value(r);
    obj_list = [obj_list,z_l];
%     z_l
    
%     disp("step 1 complete")
    
    %% step 2
    
    optimal0 = Dual_sp(cons_points,data,epsilon0,K,kappa_hat,x_hat);
    runtime_sub = runtime_sub + optimal0.runtime;
    runtime = runtime + optimal0.runtime;
    pi_hat = optimal0.pi;
    z_hat = optimal0.obj;
    z_d = optimal0.z;
    
    tic
    if z_hat <= z_u
        z_u = z_hat;
        x_opt = x_hat;
    end
    
    %% step 3
    if z_u - z_l <= error*abs(z_l)
        disp("find the optimal value")
        break
    end
    
    %% step 4
    for k = 1:K
        if r_hat(:,:,k) < z_d(k) < Inf
            count(k) = count(k) + 1;
            
            for l = 1:ll
                U(:,:,l,count(k),k) = pi_hat{l,k};
            end
            U(:,:,ll+1,count(k),k) = pi_hat{ll+1,k};
        else
            disp("The dual problem is unbounded")
        end      
    end
    
    tt = toc;
    
    iter = iter +1;
    text = ['finish ',num2str(iter),'th iteration'];
    disp(text)    
end

optimal.obj = obj_list;
optimal.q = value(x_opt);
optimal.t = toc(timerVal);
optimal.runtime = runtime;
optimal.iter = iter;
optimal.runtime_sub = runtime_sub/iter;

end