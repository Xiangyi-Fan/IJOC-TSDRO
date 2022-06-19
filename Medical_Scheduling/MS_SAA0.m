function [x,obj] =  MS_SAA0(data)


yalmip clear;

%% Parameters

N       = 8; % number of stocks
Npoints = size(data,2); % number of data points

l_l = 30; l_u = 60; % l_l <= l <= l_u
c = 100;
T = (l_l+l_u)/2*N;

% Matrix D
D = [];
for i = 1:N
    vec = zeros(N+1,1);
    vec(i) = -1;
    vec(i+1) = 1;
    D = [D;vec'];
end



for i = 1:Npoints; y{i} = sdpvar(N+1,1);end
x = sdpvar(N,1);

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = sum(x) <= T;


for i = 1:Npoints
    constraints{end+1} = D*y{i} + x >= [eye(N),zeros(N,N+1)]*data(:,i);
    constraints{end+1} = y{i} >= 0;
end


value = 0;
for i = 1:Npoints
    value = value + data(:,i)'*[zeros(N),eye(N),zeros(N,1)]'*[eye(N),zeros(N,1)]*y{i} + c*[zeros(N,1);1]'*y{i};
end
obj = value/Npoints;

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');
out = optimize([constraints{:}],obj,options);

double(obj)
double(x)

end