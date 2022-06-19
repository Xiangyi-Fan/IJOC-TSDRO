function [x,obj,out] =  pwstatic(Sp,tp,hc,budget,mu,Sigma,N,M,eps)

W = [zeros(1,2*N), 1;          %tau >= 0
     eye(2*N), zeros(2*N,1);   % y1 >= x - demand; y2 >= demand - x
     eye(2*N), zeros(2*N,1)];  % y1,y2 >= 0

x = sdpvar(N,1);
kappa = sdpvar; % beta of the cvar

hx = [0; x; -x; zeros(2*N,1)];
Tx = [zeros(1,N); -eye(N); eye(N); zeros(2*N,N)];
Tx = [zeros(size(Tx,1),N), Tx];

y1 = sdpvar(N,M,'full');
y2 = sdpvar(N,M,'full');
tau = sdpvar(M,1);

alpha = sdpvar;
beta = sdpvar(2*N,1);
Gamma = sdpvar(2*N);

K = length(hx);
rho = cell(M,1);
rho2 = cell(M,1);
rho3 = cell(M,1);

cons = [x >= 0; sum(x) <= budget];

for j = 1:M
    
    Lj = length(tp{j});
    
    %objective
    rho2{j} = sdpvar(Lj,1);
    MM = [         Gamma,          .5*(beta+Sp{j}'*rho2{j});
          .5*(beta+Sp{j}'*rho2{j})', alpha - tau(j) - tp{j}'*rho2{j}];
    cons = [cons; MM >= 0; rho2{j} >= 0];
    
    rho{j} = sdpvar(Lj,K,'full');
    cons = [cons; W*[y1(:,j);y2(:,j);tau(j)] - hx >= rho{j}'*tp{j};
                  rho{j}(:) >= 0;
                  Sp{j}'*rho{j} == Tx'];

    rho3{j} = sdpvar(Lj,1);
    cons = [cons; tau(j) + kappa - hc'*y1(:,j) >= rho3{j}'*tp{j};
                  rho3{j} >= 0;
                  Sp{j}'*rho3{j} == [y2(:,j); zeros(N,1)]];
end

obj = kappa + 1/eps * (alpha + beta'*mu + trace(Sigma*Gamma));

options = sdpsettings('verbose', 0, 'solver', 'mosek');
out = optimize(cons,obj,options);
obj = double(obj);
x = double(x);

end