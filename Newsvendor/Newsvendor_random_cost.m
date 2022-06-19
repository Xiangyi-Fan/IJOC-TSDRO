function [x,obj] =  Newsvendor_random_cost(epsilon)

yalmip clear;

% Parameters

N       = 3; % number of stocks
Npoints = 20; % number of data points

xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u
s_l = 20 ; s_u = 40; % s_l <= s <= s_u 
delta = 0.1; % risk attitude
B = 30; % total supply
epsilon1 = 0.05; % mean & convariance ambiguity set 

R = 0.1; % sqrt(N*xi_u*xi_u + N*s_u*s_u +1)
delta2 = 0.1;
epsilon2 = R^2/sqrt(Npoints)*(2+sqrt(2*log(1/delta2))); % Frobenius norm ambiguity set

b0 = max(xi_u,s_u);
b0 = max(1,b0);
a0 = min(xi_l,s_l);
a0 = min(1,a0);
delta3 = 0.1;
epsilon3 = 0.01*(b0-a0)*sqrt((log((2*N)^2/2+2*N)-log(delta3))/(2*Npoints)); % Infinite norm ambiguity set  

ee = ones(N,1);


g = [1;3;5];
P = [eye(N), zeros(N), -1*xi_l*ee;
     -1*eye(N), zeros(N), xi_u*ee;
     zeros(N), eye(N), -1*s_l*ee;
     zeros(N), -1*eye(N), s_u*ee];

data1 = NaN;
while isnan(data1)
    pd1 = makedist('Uniform','Lower',0,'Upper',2);
    mu = random(pd1,N,1);
    sigma = [1/4;1/4;1/4];
    corrmat = gallery('randcorr',3);
    data1 = MvLogNRand(mu,sigma,Npoints,corrmat);
end
for i = 1:Npoints
    for j = 1:N
        if data1(i,j) > xi_u
            data1(i,j) = xi_u;
        elseif data1(i,j) < xi_l
            data1(i,j) = xi_l;
        end
    end
end      
xi = data1';


data2 = NaN;
while isnan(data2)
    pd2 = makedist('Uniform','Lower',3.2,'Upper',3.3);
    mu = random(pd2,N,1);
    sigma = [1/4;1/4;1/4];
    corrmat = gallery('randcorr',3);
    data2 = MvLogNRand(mu,sigma,Npoints,corrmat);
end
for i = 1:Npoints
    for j = 1:N
        if data2(i,j) > s_u
            data2(i,j) = s_u;
        elseif data2(i,j) < s_l
            data2(i,j) = s_l;
        end
    end
end      
s = data2';

data = [xi;s;ones(1,Npoints)];
Sigma0 = data*data'/Npoints;
Sigma_l = (1-epsilon1)*Sigma0;
Sigma_u = (1+epsilon1)*Sigma0;

                      
S_supp = [-1*eye(N), zeros(N); eye(N), zeros(N); zeros(N), -1*eye(N); zeros(N),eye(N)];
t_supp = [-1*xi_l*ee; xi_u*ee; -1*s_l*ee; s_u*ee];
[A_supp,b_supp,~,out_supp] = polyhedron_copositive(S_supp,t_supp);

R = [b_supp'*A_supp, 0.5*b_supp'*b_supp;
     A_supp, zeros(2*N,1);
     -b_supp'*A_supp, 0.5*(2-b_supp'*b_supp)];

[n1,n2] =size(R);
S = R'*[zeros(n1-1,1);1]*[zeros(n1-1,1);1]'*R;
for i = 1:n1-1
    vec=zeros(n1,1);
    position=i;
    vec(position)=1;
    S = S-R'*vec*vec'*R;
end

%% Decision Variables

x = sdpvar(N,1);
kappa = sdpvar(1,1);
alpha = sdpvar(1,1);
Gamma_u = sdpvar(2*N+1,2*N+1,'sym');
Gamma_l = sdpvar(2*N+1,2*N+1,'sym');
gamma = sdpvar(1,1);
Q = sdpvar(2*N+1,2*N+1,'full');
theta = sdpvar(1,1);
eta = sdpvar(1,1);
Y1 = sdpvar(N,2*N+1,'full');
Y2 = sdpvar(N,2*N+1,'full');
beta = sdpvar(1,1);
V1 = sdpvar(2*N+1,2*N+1,'full');
V2 = sdpvar(2*N+1,2*N+1,'full');
V3 = sdpvar(2*N+1,2*N+1,'full');
W1 = sdpvar(2*N+1,2*N+1);
W2 = sdpvar(2*N+1,2*N+1);
W3 = sdpvar(2*N+1,2*N+1);
Sigma1 = sdpvar(4*N,4*N);
Sigma2 = sdpvar(4*N,4*N);
Sigma3 = sdpvar(4*N,4*N);
Psi1 = sdpvar(2*N+1,2*N+1);
Psi2 = sdpvar(2*N+1,2*N+1);
Psi3 = sdpvar(2*N+1,2*N+1);
Phi1 = sdpvar(4*N,n1);
Phi2 = sdpvar(4*N,n1);
Phi3 = sdpvar(4*N,n1);
tau1 = sdpvar(1,1);
tau2 = sdpvar(1,1);
tau3 = sdpvar(1,1);
tau_u = sdpvar(2*N,1);
tau_l = sdpvar(2*N,1);
n = sdpvar(N,1);
b = sdpvar(N,1);
d = sdpvar(N,1);
h = sdpvar(N,1);
M = sdpvar(4*N,N);
A = sdpvar(4*N,N);
C = sdpvar(4*N,N);
F = sdpvar(4*N,N);
lambda = sdpvar(1,1);

%% Constraints

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = sum(x) <= B;


%% 

constraints{end+1} = theta+kappa >= 0; 
constraints{end+1} = 0.5*((Q-[zeros(2*N,1);1]*g'*Y1-[zeros(N);eye(N);zeros(1,N)]*Y2)+(Q-[zeros(2*N,1);1]*g'*Y1-[zeros(N);eye(N);zeros(1,N)]*Y2)')-theta*[zeros(2*N,1);1]*[zeros(2*N,1);1]' == V2;

%% 

constraints{end+1} = beta >= 0;
constraints{end+1} = 0.5*(Q+Q')-beta*[zeros(2*N,1);1]*[zeros(2*N,1);1]' == V3;

constraints{end+1} = V1 == W1+P'*Sigma1*P;
constraints{end+1} = Sigma1(:) >= 0;
constraints{end+1} = W1 >= 0;
constraints{end+1} = V2 == W2+P'*Sigma2*P;
constraints{end+1} = Sigma2(:) >= 0;
constraints{end+1} = W2 >= 0;
constraints{end+1} = V3 == W3+P'*Sigma3*P;
constraints{end+1} = Sigma3(:)>=0;
constraints{end+1} = W3 >= 0;

%% 

constraints{end+1} = M >= 0;
constraints{end+1} = n-x >= 0;
constraints{end+1} = Y1+[eye(N),zeros(N,N+1)]-M'*P-n*[zeros(2*N,1);1]' == 0;

%%
constraints{end+1} = A >= 0;
constraints{end+1} = b+x >= 0;
constraints{end+1} = Y2-[eye(N),zeros(N,N+1)]-A'*P-b*[zeros(2*N,1);1]' == 0;

%%

constraints{end+1} = C >= 0;
constraints{end+1} = d >= 0;
constraints{end+1} = Y1-C'*P-d*[zeros(2*N,1);1]' == 0;

%%

constraints{end+1} = F>=0;
constraints{end+1} = h>=0;
constraints{end+1} = Y2-F'*P-h*[zeros(2*N,1);1]' == 0;

%% objective: min

obj = kappa+1/delta*(epsilon*norm(Q-lambda*[zeros(2*N,1);1]*[zeros(2*N,1);1]','fro') + trace(Q*Sigma0)); % Frobenius-norm + last component

%% solving and post-processing

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');

out = optimize([constraints{:}],obj,options);


end