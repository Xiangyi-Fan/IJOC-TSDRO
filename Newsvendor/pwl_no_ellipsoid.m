function [x,obj,out] =  pwl_no_ellipsoid(Sp,tp,hc,budget,mu,Sigma,N,M,eps)

W = [zeros(1,2*N), 1;          %tau >= 0
     eye(2*N), zeros(2*N,1);   % y1 >= x - demand; y2 >= demand - x
     eye(2*N), zeros(2*N,1)];  % y1,y2 >= 0

x = sdpvar(N,1);
kappa = sdpvar; % beta of the cvar
hx = [0; x; -x; zeros(2*N,1)];
Tx = [zeros(1,N); -eye(N); eye(N); zeros(2*N,N)];
Tx = [zeros(size(Tx,1),N), Tx];

K = length(hx);

Y1 = cell(M,1);
y1 = cell(M,1);
Y2 = cell(M,1);
y2 = cell(M,1);
Tau = cell(M,1);
tau = sdpvar(M,1);

for j = 1:M
    Y1{j} = sdpvar(N,2*N,'full');
    y1{j} = sdpvar(N,1);
    Y2{j} = sdpvar(N,2*N,'full');
    y2{j} = sdpvar(N,1);
    Tau{j} = sdpvar(2*N,1);
end

alpha = sdpvar;
beta = sdpvar(2*N,1);
Gamma = sdpvar(2*N);

rho = cell(M,1);
rho2 = cell(M,1);
rho3 = cell(M,1);

cons = [x >= 0; sum(x) <= budget; Gamma >= 0];

for j = 1:M
    Lj = length(tp{j});
    %objective
    rho2{j} = sdpvar(Lj,1);
    MM = [             Gamma,              .5*(beta+Sp{j}'*rho2{j}-Tau{j});
          .5*(beta+Sp{j}'*rho2{j}-Tau{j})', alpha - tau(j) - tp{j}'*rho2{j}];
    cons = [cons; MM >= 0; rho2{j} >= 0];
        
    rho{j} = sdpvar(Lj,K,'full');
    cons = [cons; W*[y1{j};y2{j};tau(j)] - hx >= rho{j}'*tp{j};
                  rho{j}(:) >= 0;
                  Sp{j}'*rho{j} == Tx' - [Y1{j}; Y2{j}; Tau{j}']'*W'];

    rho3{j} = sdpvar(Lj,1);
    YY2 = [Y2{j}; zeros(N,2*N)];
    yy2 = [y2{j}; zeros(N,1)];
    MM = [   -.5*(YY2+YY2'),       .5*(Sp{j}'*rho3{j}+Tau{j}-yy2-Y1{j}'*hc);
          .5*(Sp{j}'*rho3{j}+Tau{j}-yy2-Y1{j}'*hc)', tau(j) + kappa - hc'*y1{j} - tp{j}'*rho3{j}];
    cons = [cons; MM >= 0; rho3{j} >= 0];
    
end

obj = kappa + 1/eps * (alpha + beta'*mu + trace(Sigma*Gamma));

options = sdpsettings('verbose', 0, 'solver', 'mosek');
out = optimize(cons,obj,options);

obj = double(obj);
x = double(x);

end