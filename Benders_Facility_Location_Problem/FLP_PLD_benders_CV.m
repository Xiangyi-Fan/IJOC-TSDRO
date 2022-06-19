function [x_value] = FLP_PLD_benders_CV(data,epsilon0,f,C)

yalmip clear;

timerVal = tic;
options = sdpsettings('verbose', 0, 'dualize',false, 'solver', 'mosek');

% Parameters

I = 5;
J = 5;
Npoints = size(data,2); % number of data points
K = size(data,2); % number of partitions

delta = 1; % risk attitude


u = 1500*ones(I,1);
g = 1000*ones(J,1);

% K fold CV

epsilon = cell(K,1);
for k = 1:K
    epsilon{k,1} = epsilon0/sqrt(Npoints);
end

% Define P_k
P = cell(K,1);
cons_points = data;
for k = 1:K
    P{k,1} = FLP_PLD_partitions(cons_points,k);
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


%% step 0
z_u = Inf;
x_opt = zeros(I,1);
U = zeros(J+1,J+1,ll,0,K);
V = zeros(J+1,J+1,ll,0,K);
obj_list = [];
error = 0;
count_U = zeros(1,K);
count_V = zeros(1,K);

%% begin
runtime = 0;
logical = true;
iter = 1;
while logical
    
    %% step 1
    
    % Decision Variables
    x = binvar(I,1);
    r = sdpvar(1,1,K);
    
    % Parameters
    for l = 1:J
        T_hat{l} = ee(J+1,J+1);
    end

    for l = J+1:J+I
        l1 = l - J;
        T_hat{l} = -1*ee(J+1,J+1)*u(l1)*x(l1);
    end
    for l = J+I+1:2*J+I+I*J
        l1 = l - (J+I);
        T_hat{l} = zeros(J+1,1);
    end 
    T_hat{ll-1} = zeros(J+1,1);
    T_hat{ll} = zeros(J+1,1);
        
     
    % Constraints
    constraints1 = {};
    constraints1{end+1} = x >= 0;
    for k = 1:K
        constraints1{end+1} = r(1,1,k) >= -9999;
    end
    
    
    for k = 1:K
        for j = 1:size(U,4)
            if ~isnan(U(:,:,1,j,k))
                cons = 0;
                for l = 1:ll
                    cons = cons + 0.5*vector(U(:,:,l,j,k))'*vector(ee(J+1,J+1)*T_hat{l}'+T_hat{l}*ee(J+1,J+1)');
                end
                constraints1{end+1} = cons <= r(:,:,k);
            end
        end  
        
        for j = 1:size(V,4)
            if ~isnan(V(:,:,1,j,k))
                cons = 0;
                for l = 1:ll
                    cons = cons + 0.5*vector(V(:,:,l,j,k))'*vector(ee(J+1,J+1)*T_hat{l}'+T_hat{l}*ee(J+1,J+1)');
                end
                constraints1{end+1} = cons <= 0;
            end
        end 
    end

    % objective: min
    obj = f'*x + 1/delta*sum(r);

    % solving and post-processing
    out = optimize([constraints1{:}],obj,options);
    runtime = runtime + out.solvertime;
    z_l = double(obj);
    x_hat = value(x);
%     theta_hat = value(theta);
    r_hat = value(r);
    obj_list = [obj_list,z_l];
    
%     disp("step 1 complete")
    
    %% step 2
    
    optimal0 = Dual(data,epsilon0,K,x_hat,C);
    runtime = runtime + optimal0.runtime;
    pi_hat = optimal0.pi;
    z_d_feasibility = optimal0.feasibility;
    
    tic
    if max(z_d_feasibility) == 0
        z_hat = optimal0.obj + f'*x_hat;
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
%         disp("find the optimal value")
        tt1 = toc;
        runtime = runtime + tt1;
        break
    end
    tt1 = toc;
    runtime = runtime + tt1;
    
    %% step 4
    cut_time = zeros(K,1);
    for k = 1:K
        if z_d_feasibility(k) == 0
        
        % add optimality cut
            tic
            count_U(k) = count_U(k) + 1;  
            for l = 1:ll
                U(:,:,l,count_U(k),k) = pi_hat{l,k};
            end
            tt2 = toc;
            cut_time(k) = tt2;
%             text = [num2str(num),'th optimality cut is successfully added for N =', num2str(Npoints)];
%             disp(text)
        else
        % add feasibility cut
%             text = [num2str(k), 'th part of the dual problem is unbounded'];
%             disp(text)
            cut = Dual_feasibility_cut(data,epsilon0,k,x_hat,C);
            tic
            pi_hat0  = cut.pi;
            count_V(k) = count_V(k) + 1;
            for l = 1:ll
                V(:,:,l,count_V(k),k) = pi_hat0{l};
            end
            tt3 = toc;
            cut_time(k) =  tt3 + cut.runtime;
        end
    end      
    
    runtime = runtime + max(cut_time);
%     text  = ['finish ',num2str(iter),'th iteration'];
%     disp(text)
    iter = iter +1;

    if iter >= 30
        break
    end
    
end

x_value = value(x_opt);


end