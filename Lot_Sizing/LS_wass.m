function [x,obj,runtime] =  LS_wass(data0,epsilon0_wass)

yalmip clear;

%% Parameters

N       = 5; % number of locations
Npoints = size(data0,2); % number of data points

v_l = 40; v_u = 50; % v_l <= v <= v_u
u_l = 20 ; u_u = 40; % u_l <= u <= u_u 
C = [40;50;60;70;80];
cap = 80;

N2 = N*N;
J = 2*(N*N+N);
K2 = N*N+N;
N1 = N2+N;
b = K2+N1+J+1;

x = sdpvar(N,1); %decision variable x

% Matrix Q
Q = [eye(N2),zeros(N2,N)];

q = zeros(N2,1);
 
% Matrix T
T = [zeros(N,N2),eye(N);
     zeros(N2,K2)];
 
% Matrix h(x)
h = [   -1*x;
     zeros(N2,1)];

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

% Matrix W
W = [ A + B;
     eye(N2)];

% Matrix S
ee1 = ones(N,1);
ee2 = ones(N2,1);
S = [eye(N2), zeros(N2,N);
     -1*eye(N2), zeros(N2,N);
     zeros(N,N2), eye(N);
     zeros(N,N2), -1*eye(N)];

tt = [v_u*ee2;-1*v_l*ee2;u_u*ee1;-1*u_l*ee1];


% Modified Matrix
QQ = [Q;S];
qq = [q;-1*tt];
TT = [T;zeros(J,K2)];
hh = [h;zeros(J,1)];
WW = [W,zeros(N1,J);
      zeros(J,N2),-1*eye(J)];

data0 = data0(1:(N2+N),:);

%% Decision Variables

lambda = sdpvar(1,1);
for i = 1:Npoints; s{i} = sdpvar(1);end
for i = 1:Npoints; phi{i} = sdpvar(N2+J,1);end
for i = 1:Npoints; psi{i} = sdpvar(N2+J,1);end
for i = 1:Npoints; V0{i} = sdpvar(b,b);end
for i = 1:Npoints; P0{i}= sdpvar(b,b);end
for i = 1:Npoints; W0{i} = sdpvar(b,b);end

%% Constraints

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = x <= cap;
constraints{end+1} = lambda >= 0;

for i = 1:Npoints
        %% 
        constraints{end+1} = [lambda*eye(K2)+QQ'*diag(phi{i})*QQ,     -0.5*TT'-QQ'*diag(phi{i})*WW', -1*lambda*data0(:,i)-0.5*QQ'*psi{i};
                              -0.5*TT-WW*diag(phi{i})*QQ,              WW*diag(phi{i})*WW',                  0.5*(WW*psi{i}-hh);
                              (-1*lambda*data0(:,i)-0.5*QQ'*psi{i})',  0.5*(WW*psi{i}-hh)',                         s{i}]                 == V0{i};
               
        
        %% SDP Approximation
        
        constraints{end+1} = V0{i} == P0{i} + W0{i};
        constraints{end+1} = P0{i} >= 0;
        constraints{end+1} = W0{i}(:) >= 0;
        
end

%% objective: min
sum_s = 0;
for i = 1:Npoints
    sum_sub = 0;
    for j = 1:(N2+J)
        sum_sub = sum_sub + phi{i}(j)*qq(j)*qq(j);
    end
    sum_s = sum_s + s{i} + qq'*psi{i} - lambda*(data0(:,i)')*data0(:,i) + sum_sub;
end

obj = C'*x + epsilon0_wass^2*lambda + 1/Npoints*sum_s;

%% solving and post-processing

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');

out = optimize([constraints{:}],obj,options);
runtime = out.solvertime;

end