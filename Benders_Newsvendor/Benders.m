function [optimal] = Benders(cons_points,data,epsilon0,K)

yalmip clear;

timerVal = tic;
options = sdpsettings('verbose', 0, 'solver', 'mosek');

%% Parameters

N       = 5; % number of stocks
Npoints = size(data,2); % number of data points

xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u
s_l = 0; s_u = 50; % s_l <= s <= s_u 
delta = 0.1; % risk attitude
B = 30; % total supply
g = [5;6;7;8;9];

% Define epsilon_k
epsilon = cell(K,1);
for k = 1:K
    epsilon{k,1} = epsilon0/sqrt(Npoints);
end

% Define P_k
P = cell(K,1);
for k = 1:K
    P{k,1} = PLD_partitions(cons_points,k,K);
end

vector = @(x) x(:); %transfer a matrix to a vector

%% Index
aa = 2*N+1;
bb = N;
cc = 4*N;
ll = 2;
nn = 2;
mm = 1;

%% step 0
z_u = Inf;
x_opt = zeros(N,1);
U1 = zeros(aa,aa,ll,0,K);
U2 = zeros(cc,aa,1,0,K);
obj_list = [];
error = 0.05;
count = zeros(1,K);

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
    T = cell(ll,1);
    T{1} = [zeros(2*N,1);-1*kappa];
    T{2} = zeros(aa,1);
    G = [-1*eye(bb),zeros(bb),x;
         eye(bb),zeros(bb),-1*x;
         zeros(bb,aa);
         zeros(bb,aa)];
       
     
    % Constraints
    constraints1 = {};
    
    constraints1{end+1} = kappa >= -9999;
    for k = 1:K
        constraints1{end+1} = r(1,1,k) >= -9999;
    end
    
    constraints1{end+1} = x >= 0;
    constraints1{end+1} = sum(x) <= B;
    
    for k = 1:K
        for j = 1:size(U1,4)
            if ~isnan(U1(:,:,1,j,k))
                cons = -1*trace(G'*U2(:,:,1,j,k));
                for l = 1:ll
                    cons = cons + 0.5*trace(U1(:,:,l,j,k)'*(ee*T{l}'+T{l}*ee'));
                end
                constraints1{end+1} = cons <= r(:,:,k);
            end
        end  
    end

    % objective: min
    obj = kappa;
    for k = 1:K
        I_k =  News_PLD_nums(P{k,1}, data);
        obj = obj + 1/delta*r(:,:,k)*I_k/Npoints;
    end

    % solving and post-processing
    out = optimize([constraints1{:}],obj,options);
    runtime = runtime + out.solvertime;
    z_l = double(obj);
    x_hat = value(x);
    kappa_hat = value(kappa);
    r_hat = value(r);
    obj_list = [obj_list,z_l];
    
    %% step 2
    
    optimal0 = Dual_sp(cons_points,data,epsilon0,K,kappa_hat,x_hat);
    runtime_sub = runtime_sub + optimal0.runtime;
    runtime = runtime + optimal0.runtime;
    pi_hat = optimal0.pi;
    z_hat = optimal0.obj;
    z_d = optimal0.z;
    z_d_feasibility = optimal0.feasibility;
    
    tic
    if z_hat <= z_u
        z_u = z_hat;
        x_opt = x_hat;
    end
    
    %% step 3
    if z_u - z_l <= error*abs(z_l)
        tt = toc;
        disp("find the optimal value")
        break
    end
    
    %% step 4
    for k = 1:K
        if z_d_feasibility(k) == 0
            count(k) = count(k) + 1;
            
            for l = 1:ll
                U1(:,:,l,count(k),k) = pi_hat{l,k};
            end
            U2(:,:,1,count(k),k) = pi_hat{ll+1,k};
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