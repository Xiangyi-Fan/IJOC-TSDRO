function [obj,Y1,Y2] = Newsvendor_PLD_cop_fix(data)

yalmip clear;

%% Parameters

N       = 2; % number of stocks
Npoints = size(data,2); % number of data points
K = size(data,2); % number of partitions

xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u
delta = 0.1; % risk attitude
B = 20; % total supply
g = [2;4];
s = [10;12];

data = [data;ones(1,Npoints)];

% Define P_k
P = cell(K,1);
cons_points = data;
for k = 1:K
    P{k,1} = PLD_partitions_fix(cons_points,k);
end

u0 = cell(K,1);
for k = 1:K
    u0{k,1} = zeros(size(data,1),1);
    for i = 1:Npoints
        if P{k,1}*data(:,i) >= 0
            u0{k,1} = u0{k,1} + data(:,i);
        end
    end
    u0{k,1} = u0{k,1}/Npoints;
end

x = [3.96;6.25];
%% Decision Variables

l = size(P{1},1);
% x = sdpvar(N,1);
kappa = sdpvar(1,1);
for k = 1:K; y{k} = sdpvar(N+1,1);end
for k = 1:K; Y1{k} = sdpvar(N,N+1);end
for k = 1:K; Y2{k} = sdpvar(N,N+1);end
for k = 1:K; alpha{k} = sdpvar(l,1);end
for k = 1:K; beta{k} = sdpvar(1,1);end
for k = 1:K; gamma{k} = sdpvar(l,1);end
for k = 1:K; eta{k} = sdpvar(1,1);end
for k = 1:K; n{k} = sdpvar(N,1);end
for k = 1:K; b{k} = sdpvar(N,1);end
for k = 1:K; d{k} = sdpvar(N,1);end
for k = 1:K; h{k} = sdpvar(N,1);end
for k = 1:K; M{k} = sdpvar(l,N);end
for k = 1:K; A{k} = sdpvar(l,N);end
for k = 1:K; C{k} = sdpvar(l,N);end
for k = 1:K; F{k} = sdpvar(l,N);end

%% Constraints

constraints = {};

for k = 1:K
    %% 

    constraints{end+1} = M{k} >= 0;
    constraints{end+1} = n{k}-x >= 0;
    constraints{end+1} = Y1{k}+[eye(N),zeros(N,1)]-M{k}'*P{k}-n{k}*[zeros(N,1);1]' == 0;

    %%
    
    constraints{end+1} = A{k} >= 0;
    constraints{end+1} = b{k}+x >= 0;
    constraints{end+1} = Y2{k}-[eye(N),zeros(N,1)]-A{k}'*P{k}-b{k}*[zeros(N,1);1]' == 0;

    %% 
    
    constraints{end+1} = alpha{k} >= 0;
    constraints{end+1} = beta{k}+kappa >= 0; 
    constraints{end+1} = y{k} - Y1{k}'*g - Y2{k}'*s -  P{k}'*alpha{k} - beta{k}*[zeros(N,1);1] == 0;

    %%

    constraints{end+1} = C{k} >= 0;
    constraints{end+1} = d{k} >= 0;
    constraints{end+1} = Y1{k}-C{k}'*P{k}-d{k}*[zeros(N,1);1]' == 0;

    %%

    constraints{end+1} = F{k} >= 0;
    constraints{end+1} = h{k} >= 0;
    constraints{end+1} = Y2{k}-F{k}'*P{k}-h{k}*[zeros(N,1);1]' == 0;

    %% 
    constraints{end+1} = gamma{k} >= 0;
    constraints{end+1} = eta{k} >= 0; 
    constraints{end+1} = y{k} - P{k}'*gamma{k} - eta{k}*[zeros(N,1);1] == 0;
    
end

%% objective: min

value = cell(K,1);
values = 0;
for k = 1:K
    value{k} = y{k}'*u0{k};
    values = values + value{k};
end
     
obj = kappa+1/delta*values; % Infinite-norm + last component

%% solving and post-processing

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');

out = optimize([constraints{:}],obj,options);

end