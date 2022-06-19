function [optimal] = Benders_general(cons_points,data,epsilon0,K)

yalmip clear;

timerVal = tic;
options = sdpsettings('verbose', 0, 'dualize',false, 'solver', 'mosek');

% Parameters

N       = 5; % number of locations
Npoints = size(data,2); % number of data points

v_l = 40; v_u = 50; % v_l <= v <= v_u
u_l = 20 ; u_u = 40; % u_l <= u <= u_u 
C = [40;50;60;70;80];
cap = 80;


epsilon = cell(K,1);
for k = 1:K
    epsilon{k,1} = epsilon0/sqrt(Npoints);
end

pho2 = 0.1;
epsilon_p = 1/Npoints*(K - 1 + 2*sqrt((1-K)*log(pho2)) - 2*log(pho2));


% Define P_k
P = cell(K,1);
for k = 1:K
    P{k,1} = LS_PLD_partitions(cons_points,k);
end

Omega = cell(K,1);
for k = 1:K
    I_k =  LS_PLD_nums(P{k,1}, data);
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

% Matrix A and B
A = [];
B = [];
for i = 1:N
    vec=zeros(N,1);
    positions1=i;
    vec(positions1)=1;
    A0 = repmat(vec,N,1);
    A = [A;A0'];
    B0 = zeros(N*N,1);
    positions2 = [(i-1)*N+1:i*N];
    B0(positions2)=-1;
    B = [B;B0'];
end

vector = @(x) x(:); %transfer a matrix to a vector

%% Index
aa = N*N+N+1;
bb = N*N+N;
cc = N*N;

%% step 0
z_u = Inf;
x_opt = zeros(N,1);
U = zeros(bb,aa,1,0,K);
V = zeros(bb,aa,1,0,K);
obj_list = [];
error = 0.05;
count_U = zeros(1,K);
count_V = zeros(1,K);

%% begin
runtime = 0;
runtime_sub = 0;
logical = true;
iter = 0;
while logical
    
    %% step 1
    
    % Decision Variables
    x = sdpvar(N,1);
    ss = sdpvar(1,1,K);
    eta = sdpvar(1,1);
    omega = sdpvar(1,1);
    rr = sdpvar(1,1,K);
    
    % Parameters
    G = [zeros(N,cc),eye(N),-1*x;
         zeros(cc,aa)];    
     
    % Constraints
    constraints1 = {};
    
    for k = 1:K
        constraints1{end+1} = ss(1,1,k) >= -9999;
    end
    
    constraints1{end+1} = x >= 0;
    constraints1{end+1} = x <= cap;
    
    for k = 1:K
        for j = 1:size(U,4)
            if ~isnan(U(:,:,1,j,k))
                cons = -1*trace(G'*U(:,:,1,j,k));
                constraints1{end+1} = cons <= ss(:,:,k);
            end
        end  
        for j = 1:size(V,4)
            if ~isnan(V(:,:,1,j,k))
                cons = -1*trace(G'*V(:,:,1,j,k));
                constraints1{end+1} = cons <= 0;
            end
        end 
    end
    
    
    for k = 1:K
        constraints1{end+1} = ss(1,1,k) + eta <= omega;
        constraints1{end+1} = sqrt(4*rr(1,1,k)^2 + (ss(1,1,k)+eta)^2) <= 2*omega - ss(1,1,k) - eta;
    end
    
    constraints1{end+1} = omega >= 0;

    % objective: min
  
    obj = C'*x + epsilon_p*omega - eta + 2*omega*1; % Infinite-norm + last component
    for k = 1:K
        I_k =  LS_PLD_nums(P{k,1}, data);
        obj = obj - 2*rr(1,1,k)*I_k/Npoints;   
    end

    % solving and post-processing
    out = optimize([constraints1{:}],obj,options);
    runtime = runtime + out.solvertime;
    z_l = double(obj);
    x_hat = value(x);
    obj_list = [obj_list,z_l];
    
%     disp("step 1 complete")
    
    %% step 2
    
    optimal0 = Dual(cons_points,data,epsilon0,K,x_hat);
    runtime_sub = runtime_sub + optimal0.runtime;
    runtime = runtime + optimal0.runtime;
    pi_hat = optimal0.pi;
    z_d_feasibility = optimal0.feasibility;
    
    
    ss = sdpvar(1,1,K);
    eta = sdpvar(1,1);
    omega = sdpvar(1,1);
    rr = sdpvar(1,1,K);
    
    % Parameters
    G_hat = [zeros(N,cc),eye(N),-1*x_hat;
             zeros(cc,aa)];
     
    % Constraints
    constraints2 = {};
    
    for k = 1:K
        constraints2{end+1} = -1*vector(pi_hat{1,k})'*vector(G_hat) <= ss(1,1,k);
        constraints2{end+1} = ss(1,1,k) + eta <= omega;
        constraints2{end+1} = sqrt(4*rr(1,1,k)^2 + (ss(1,1,k)+eta)^2) <= 2*omega - ss(1,1,k) - eta;
    end
    
    constraints2{end+1} = omega >= 0;
    
    % objective: min
  
    obj2 = C'*x_hat + epsilon_p*omega - eta + 2*omega*1; % Infinite-norm + last component
    for k = 1:K
        I_k =  LS_PLD_nums(P{k,1}, data);
        obj2 = obj2 - 2*rr(1,1,k)*I_k/Npoints;   
    end

    % solving and post-processing
    options2 = sdpsettings('verbose', 0, 'dualize',true, 'solver', 'mosek');
    out2 = optimize([constraints2{:}],obj2,options2);
    runtime = runtime + out2.solvertime;
    
    tic
    if max(z_d_feasibility) == 0
        z_hat = double(obj2);
    else
        z_hat = Inf;
    end
           
    if z_hat <= z_u
        z_u = z_hat;
        x_opt = x_hat;
    end
%     disp("step 2 complete")


    %% step 3
    if z_u - z_l <= error*abs(z_l)
        disp("find the optimal value")
        tt1 = toc;
        break
    end
    tt1 = toc;
    
    %% step 4
    cut_time = zeros(K,1);
    for k = 1:K
        if z_d_feasibility(k) == 0
        
        % add optimality cut
            tic
            count_U(k) = count_U(k) + 1;            
            U(:,:,1,count_U(k),k) = pi_hat{1,k};
            tt2 = toc;
            cut_time(k) = 0;
        else
        % add feasibility cut
            cut = Dual_feasibility_cut(cons_points,data,epsilon0,k,x_hat);
            tic
            pi_hat0  = cut.pi;
            count_V(k) = count_V(k) + 1;
            V(:,:,1,count_V(k),k) = pi_hat0;
            tt3 = toc;
            cut_time(k) =  cut.runtime;
        end
    end      
    
    runtime = runtime + max(cut_time);
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