function [R_k] =  News_PLD_radius(P_k, cons_data, k)
yalmip clear;

N = 5; % number of locations
b = 2*N + 1;
data_k = cons_data(:,k);

%% Decision variables
xi = sdpvar(b,1);

%% Constraints
constraints = {};
constraints{end+1} = P_k*xi >= 0;

%% Objective
obj = norm(xi - data_k, 1);
obj = -1*obj;

%% solving and post-processing

options = sdpsettings('dualize', 0, 'verbose', 0, 'solver', 'mosek');

out = optimize([constraints{:}],obj,options);
runtime = out.solvertime;
R_k = norm(double(xi) - data_k, 2);

end