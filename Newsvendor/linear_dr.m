function [x,obj,out,kappa] =  linear_dr(S_supp,t_supp,hc,budget,mu,Sigma,N,A,b,eps)

yalmip clear;

Sp = S_supp;
tp = t_supp;

W = [zeros(1,2*N), 1;          %tau >= 0
     eye(2*N), zeros(2*N,1);   % y1 >= x - demand; y2 >= demand - x
     eye(2*N), zeros(2*N,1)];  % y1,y2 >= 0

x = sdpvar(N,1);
kappa = sdpvar; % beta of the cvar

hx = [0; x; -x; zeros(2*N,1)];
Tx = [zeros(1,N); -eye(N); eye(N); zeros(2*N,N)];
Tx = [zeros(size(Tx,1),N), Tx];

K = length(hx);

Y1 = sdpvar(N,2*N,'full');
y1 = sdpvar(N,1);
Y2 = sdpvar(N,2*N,'full');
y2 = sdpvar(N,1);
Tau = sdpvar(2*N,1);
tau = sdpvar;

alpha = sdpvar;
beta = sdpvar(2*N,1);
Gamma = sdpvar(2*N);

lambda = sdpvar;
lambda2 = sdpvar;

cons = [x >= 0; sum(x) <= budget; Gamma >= 0; lambda >= 0; lambda2 >= 0];

Lj = length(tp);

%objective
rho2 = sdpvar(Lj,1);
JJ = [A^2, A'*b; b'*A, b'*b-1 ];
MM = [      Gamma,              .5*(beta+Sp'*rho2-Tau);
      .5*(beta+Sp'*rho2-Tau)', alpha - tau - tp'*rho2];
cons = [cons; MM + lambda*JJ >= 0; rho2 >= 0];

rho = sdpvar(Lj,K,'full');
cons = [cons; W*[y1;y2;tau] - hx >= rho'*tp;
              rho(:) >= 0;
              Sp'*rho == Tx' - [Y1; Y2; Tau']'*W'];

rho3 = sdpvar(Lj,1);
YY2 = [Y2; zeros(N,2*N)];
yy2 = [y2; zeros(N,1)];
MM = [   -.5*(YY2+YY2'),       .5*(Sp'*rho3+Tau-yy2-Y1'*hc);
      .5*(Sp'*rho3+Tau-yy2-Y1'*hc)', tau + kappa - hc'*y1 - tp'*rho3];
cons = [cons; MM + lambda2*JJ >= 0; rho3 >= 0];
    

obj = kappa + 1/eps * (alpha + beta'*mu + trace(Sigma*Gamma));

options = sdpsettings('verbose', 0, 'solver', 'mosek');
out = optimize(cons,obj,options);

obj = double(obj);
x = double(x);
kappa = double(kappa);

end