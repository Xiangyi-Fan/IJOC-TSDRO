function [x,obj] =  SAA(data)

yalmip clear;

% Parameters

N       = 3; % number of stocks
Npoints = size(data,2); % number of data points

xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u
s_l = 20; s_u = 40; % s_l <= s <= s_u 
delta = 0.1; % risk attitude
B = 30; % total supply

g = [1;3;5];
  
xi = data(1:N,:);
s = data(N+1:2*N,:);

%% Decision Variables

kappa = sdpvar(1,1);
x = sdpvar(N,1);
tau = sdpvar(1,1,Npoints);
y1 = sdpvar(N,1,Npoints);
y2 = sdpvar(N,1,Npoints);
% alpha = sdpvar(1,1);

%% Constraints

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = sum(x) <= B; 

%% 

for i = 1:Npoints
    % constraints{end+1} = alpha >= tau(:,:,i);
    constraints{end+1} = tau(:,:,i) >= 0;
    constraints{end+1} = tau(:,:,i) >= g'*y1(:,:,i)+s(:,i)'*y2(:,:,i)-kappa;
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