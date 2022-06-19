function [x,obj] =  LS_SAA0(data)


yalmip clear;

%% Parameters

N       = 5; % number of stocks
Spoints = size(data,2); % number of data points
a = N*N;

C = [50;60;70;80;90]*0.8;
cap = 80;

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

T = [zeros(N,a),eye(N),zeros(N,1)];
R = [eye(a);zeros(N+1,a)];

for i = 1:Spoints; y{i} = sdpvar(a,1);end
x = sdpvar(N,1);

constraints = {};
constraints{end+1} = x >= 0;
constraints{end+1} = x <= cap;

for i = 1:Spoints
    constraints{end+1} = A*y{i} + B*y{i} + x >= T*data(:,i);
    constraints{end+1} = y{i} >= 0;
end

obj = C'*x;
value = 0;
for i = 1:Spoints
    value = value + data(:,i)'*R*y{i};
end
obj = obj + value/Spoints;

options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'mosek');
out = optimize([constraints{:}],obj,options);

double(obj)
double(x)

end

