function [x,obj,y1,y2] =  SAA_fix2(data)

yalmip clear;

% Parameters

N       = 2; % number of stocks
Npoints = size(data,2); % number of data points


delta = 0.1; % risk attitude
B = 20; % total supply
g = [2;4];
xi = data;
s = [10;12];

%% Decision Variables

kappa = sdpvar(1,1);
x = sdpvar(N,1);
tau = sdpvar(1,1,Npoints);
y1 = sdpvar(N,1,Npoints);
y2 = sdpvar(N,1,Npoints);


%% Constraints

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = sum(x) <= B; 

%% 

for i = 1:Npoints
    constraints{end+1} = tau(:,:,i) >= 0;
    constraints{end+1} = tau(:,:,i) >= g'*y1(:,:,i)+s'*y2(:,:,i)-kappa;
    constraints{end+1} = y1(:,:,i) >= 0;
    constraints{end+1} = y1(:,:,i) >= x - xi(:,i);
    constraints{end+1} = y2(:,:,i) >= 0;
    constraints{end+1} = y2(:,:,i) >= xi(:,i) - x;
end


%% objective: min

obj = kappa+1/delta*1/Npoints*sum(tau);

%% solving and post-processing

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');

out = optimize([constraints{:}],obj,options);


end