% Newsvendor Model
% The second stage problem is:
% Z(x,xi) = inf tau st tau >= h'*y1 + s'*y2 - beta, tau >= 0,
%                      y1 >= x-xi, y2 >= xi-x, y1, y2 >= 0;
% Demand and stockout cost are uncertain

clear all;
rng(15)

N = 7; % number of products
M = 4; % number of partitions
nexperiments = 1;

eps = 0;

imp_ps = zeros(nexperiments,1);
imp_pwl2 = zeros(nexperiments,1);
imp_pwl3 = zeros(nexperiments,1);
imp_ldr = zeros(nexperiments,1);

obj_ps = zeros(nexperiments,1);
obj_pwl = zeros(nexperiments,1);
obj_pwl2 = zeros(nexperiments,1);
obj_pwl3 = zeros(nexperiments,1);
obj_ldr = zeros(nexperiments,1);

time_ps = zeros(nexperiments,1);
time_pwl = zeros(nexperiments,1);
time_pwl2 = zeros(nexperiments,1);
time_pwl3 = zeros(nexperiments,1);
time_ldr = zeros(nexperiments,1);

xps = cell(nexperiments,1);
xpwl = cell(nexperiments,1);
xpwl2 = cell(nexperiments,1);
xpwl3 = cell(nexperiments,1);
xldr = cell(nexperiments,1);

out_ps = cell(nexperiments,1);
out_pwl = cell(nexperiments,1);
out_pwl2 = cell(nexperiments,1);
out_pwl3 = cell(nexperiments,1);
out_ldr = cell(nexperiments,1);

for i = 1:nexperiments
    
    budget = 30;
    hc = 2*rand(N,1); % holding cost
    

    xil = 0;  % demand lb,ub
    xiu = 10;
    muxi = 2*rand(N,1);
    
    scl = 8;% stockout cost lb,ub
    scu = 12;
    musc = (scl+scu)/2 * ones(N,1);

    mu = [musc; muxi];
    C = gallery('randcorr',2*N);
    
    sigma = [ones(N,1)/2; muxi/4]; 
    Sigma = diag(sigma)*diag(sigma) + mu*mu';

    S_supp = [-eye(2*N); eye(2*N)];
    t_supp = [-scl*ones(N,1);-xil*ones(N,1); scu*ones(N,1); xiu*ones(N,1)];
    points = [scl + (scu-scl)*rand(M,N), xil + (xiu-xil)*rand(M,N)];
    
    [A_supp,b_supp,~,out_supp] = polyhedron_copositive(S_supp,t_supp);
    
    [Sp,tp,A,b,ellip_time] = get_ellipsoids(S_supp,t_supp,points);
    

    [xps{i},obj_ps(i),out_ps{i}] =  pwstatic(Sp,tp,hc,budget,mu,Sigma,N,M,eps);
    [xpwl{i},obj_pwl(i),out_pwl{i}] =  pw_linear(Sp,tp,hc,budget,mu,Sigma,N,M,A,b,eps);
    [xpwl2{i},obj_pwl2(i),out_pwl2{i}] =  pwl_double_rad(Sp,tp,hc,budget,mu,Sigma,N,M,A,b,eps);
    [xldr{i},obj_ldr(i),out_ldr{i}] =  linear_dr(S_supp,t_supp,hc,budget,mu,Sigma,N,A_supp,b_supp,eps);

    imp_ps(i) = 100*(obj_ps(i) - obj_pwl(i))/obj_pwl(i);
    imp_pwl2(i) = 100*(obj_pwl2(i) - obj_pwl(i))/obj_pwl(i);
    imp_ldr(i) = 100*(obj_ldr(i) - obj_pwl(i))/obj_pwl(i);
    
    [i, imp_ps(i), imp_pwl2(i), imp_ldr(i)];
    [mean(imp_ps(1:i)), mean(imp_pwl2(1:i)),  mean(imp_ldr(1:i))];
    obj_ps(i)
    obj_pwl(i)
    obj_pwl2(i)
    obj_ldr(i)
 
    time_ps(i) = out_ps{i}.solvertime;
    time_pwl(i) = out_pwl{i}.solvertime + ellip_time/4;
    time_pwl2(i) = out_pwl2{i}.solvertime + ellip_time/4;
    time_ldr(i) = out_ldr{i}.solvertime + out_supp.solvertime;

end
