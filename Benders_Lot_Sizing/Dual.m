function [optimal] = Dual(cons_points,data,epsilon0,K,x_hat)

yalmip clear;

timerVal = tic;
options = sdpsettings('verbose', 0, 'dualize',true, 'solver', 'mosek');

% Parameters

N       = 5; % number of locations
Npoints = size(data,2); % number of data points


epsilon = cell(K,1);
for k = 1:K
    epsilon{k,1} = epsilon0/sqrt(Npoints);
end

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


% Index
aa = N*N+N+1;
bb = N*N+N;
cc = N*N;


% Parameters
D = [eye(cc),zeros(cc,N+1)];
F = [A+B; eye(cc)];
G_hat = [zeros(N,cc),eye(N),-1*x_hat;
         zeros(cc,aa)];
  
     
z_d = zeros(K,1);
pi_hat = cell(1,K);
runtime = zeros(K,1);
feasibility = zeros(K,1);

for k = 1:K
    
    % Decision Variables
    M = sdpvar(bb,aa,'full');
    H = sdpvar(aa,aa,'full');
    O = sdpvar(aa,aa);
    
    % Constraints
    constraints2 = {}; 
    
    constraints2{end+1} = O >= 0;
    
    constraints2{end+1} = P{k}*O*ee(aa,aa) >= 0 ;    

    constraints2{end+1} = Omega{k}+ H - O == 0;

    constraints2{end+1} =  ee(aa,aa)'*O*ee(aa,aa) == 1;

         
    constraints2{end+1} = norm(H,'fro') <= epsilon{k};
    constraints2{end+1} = D*Omega{k}'+D*H'+F'*M == 0;
    
    matrix2 = P{k}*M';
    constraints2{end+1} = matrix2(:) <= 0;
    constraints2{end+1} = M*ee(aa,aa) <=0 ;
        
    % objective: sup
    cons3 = -1*vector(M)'*vector(G_hat);
    obj2 = cons3;
    obj2 = -1*obj2;
    
    % solving and post-processing
    out2 = optimize([constraints2{:}],obj2,options);
    runtime(k) = out2.solvertime;
    feasibility(k) = out2.problem;
    
    z_d(k) = double(-1*obj2);
    pi_hat{1,k} = value(M); 
end

z_hat = 0;
for k = 1:K
    I_k =  LS_PLD_nums(P{k,1}, data);
    z_hat = z_hat + I_k/Npoints*z_d(k);   
end 
   
optimal.obj = z_hat;
optimal.z = z_d;
optimal.pi = pi_hat;
optimal.t = toc(timerVal);
optimal.runtime = max(runtime);
optimal.feasibility = feasibility;

end