function [cut] = Dual_feasibility_cut(cons_points,data,epsilon0,k,x_hat)

yalmip clear;

options = sdpsettings('verbose', 0, 'dualize',false, 'solver', 'mosek');

% Parameters

N       = 5; % number of locations
Npoints = size(data,2); % number of data points


v_l = 40; v_u = 50; % v_l <= v <= v_u
u_l = 20 ; u_u = 40; % u_l <= u <= u_u 
C = [40;50;60;70;80];
cap = 80;

epsilon = epsilon0/sqrt(Npoints);

% Define P_k
P = LS_PLD_partitions(cons_points,k);

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
  
     
% Decision Variables
M = sdpvar(bb,aa,'full');
H = sdpvar(aa,aa,'full');
O = sdpvar(aa,aa);
    
% Constraints
constraints2 = {};  

constraints2{end+1} = O >= 0;
    
constraints2{end+1} = P*O*ee(aa,aa) >= 0 ;    

constraints2{end+1} = H - O == 0;

constraints2{end+1} =  ee(aa,aa)'*O*ee(aa,aa) == 0;
    
constraints2{end+1} = norm(H,'fro') <= epsilon;
constraints2{end+1} = D*H'+F'*M == 0;
    
matrix2 = P*M';
constraints2{end+1} = matrix2(:) <= 0;
constraints2{end+1} = M*ee(aa,aa) <=0 ;

constraints2{end+1} = H(:) >= -1*ones((aa)^2,1) ;
constraints2{end+1} = H(:) <= 1*ones((aa)^2,1) ;
constraints2{end+1} = O(:) >= -1*ones((aa)^2,1) ;
constraints2{end+1} = O(:) <= 1*ones((aa)^2,1) ;
constraints2{end+1} = M(:) >= -1*ones(bb*aa,1) ;
constraints2{end+1} = M(:) <= 1*ones(bb*aa,1) ;
        
% objective: sup
cons3 = -1*vector(M)'*vector(G_hat);
obj2 = cons3;
obj2 = -1*obj2;
    
% solving and post-processing
out2 = optimize([constraints2{:}],obj2,options);


z_d = double(-1*obj2);
pi_hat = value(M); 
   
cut.pi = pi_hat;
cut.runtime = out2.solvertime;

end