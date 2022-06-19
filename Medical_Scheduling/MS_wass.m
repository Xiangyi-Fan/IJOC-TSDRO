function [x,obj] =  MS_wass(data0,epsilon0_wass)

yalmip clear;

%% Parameters

N       = 8; % number of patients
Npoints = size(data0,2); % number of data points
 
l_l = 20; l_u = 100; % l_l <= l <= l_u
pi_l = 1 ; pi_u = 10; % pi_l <= pi <= pi_u 
c = 200;
T0 = (l_l+l_u)/2*N;
delta = 0.1; % risk attitude

x = sdpvar(N,1); %decision variable x

% Matrix Q
Q = [zeros(N),eye(N);
     zeros(1,N),zeros(1,N)];

q = [zeros(N,1);c];
 
% Matrix T
T = [eye(N),zeros(N,N);
     zeros(N+1,2*N)];
 
% Matrix h(x)
h = [-1*x;zeros(N+1,1)];

% Matrix S
ee = ones(N,1);
S = [eye(N), zeros(N);
     -1*eye(N), zeros(N);
     zeros(N), eye(N);
     zeros(N), -1*eye(N)];

tt = [l_u*ee;-1*l_l*ee;pi_u*ee;-1*pi_l*ee];

% Matrix D
D = [];
for i = 1:N
    vec = zeros(N+1,1);
    vec(i) = -1;
    vec(i+1) = 1;
    D = [D;vec'];
end

% Matrix W
W = [D;eye(N+1)];

% Modified Matrix
QQ = [Q;S];
qq = [q;-1*tt];
TT = [T;zeros(4*N,2*N)];
hh = [h;zeros(4*N,1)];
WW = [W,zeros(2*N+1,4*N);
      zeros(4*N,N+1),-1*eye(4*N)];
  
alpha = cell(2);
alpha{1} = 1/delta;
alpha{2} = 0;

data0 = data0(1:2*N,:);

%% Decision Variables
N2 = N+1;
J = 4*N;
K = 2*N;
N1 = 2*N+1;
b = K+N1+J+1;
theta = sdpvar(1,1);
lambda = sdpvar(1,1);
for i = 1:Npoints; s{i} = sdpvar(1);end
for i = 1:Npoints; for t = 1:2; kappa{i}{t} = sdpvar(1);end;end
for i = 1:Npoints; for t = 1:2; phi{i}{t} = sdpvar(N2+J,1);end;end
for i = 1:Npoints; for t = 1:2; psi{i}{t} = sdpvar(N2+J,1);end;end
for i = 1:Npoints; for t = 1:2; V0{i}{t} = sdpvar(b,b);end;end
for i = 1:Npoints; for t = 1:2; P0{i}{t} = sdpvar(b,b);end;end
for i = 1:Npoints; for t = 1:2; W0{i}{t} = sdpvar(b,b);end;end

%% Constraints

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = sum(x) <= T0;
constraints{end+1} = lambda >= 0;

for i = 1:Npoints
    for t = 1:2
        %% 
        constraints{end+1} = [lambda*eye(K)+QQ'*diag(phi{i}{t})*QQ,     -0.5*alpha{t}*TT'-QQ'*diag(phi{i}{t})*WW', -1*lambda*data0(:,i)-0.5*QQ'*psi{i}{t};
                              -0.5*alpha{t}*TT-WW*diag(phi{i}{t})*QQ,   WW*diag(phi{i}{t})*WW',                    0.5*(WW*psi{i}{t}-alpha{t}*hh);
                              (-1*lambda*data0(:,i)-0.5*QQ'*psi{i}{t})', 0.5*(WW*psi{i}{t}-alpha{t}*hh)',           s{i}+kappa{i}{t}] == V0{i}{t};
               
        %%
        sum_term = 0;
        for j = 1:N2+J
            sum_term = sum_term + phi{i}{t}(j)*(qq(j)^2);
        end
        constraints{end+1} = kappa{i}{t} == alpha{t}*theta - qq'*psi{i}{t} + lambda*norm(data0(:,i),2)*norm(data0(:,i),2) - sum_term;
        
        %% SDP Approximation
        
        constraints{end+1} = V0{i}{t} == P0{i}{t} + W0{i}{t};
        constraints{end+1} = P0{i}{t} >= 0;
        constraints{end+1} = W0{i}{t}(:) >= 0;
        
    end
end

%% objective: min
sum_s = 0;
for i = 1:Npoints
    sum_s = sum_s + s{i};
end

obj = theta + epsilon0_wass^2*lambda + 1/Npoints*sum_s;

%% solving and post-processing

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');

out = optimize([constraints{:}],obj,options);

end