yalmip clear;
A = [-1 2 0;-3 -4 1;0 0 -2];
P = sdpvar(3,3);
F = [P >= 0, A'*P+P*A <= 0, trace(P)==1];
options = sdpsettings('dualize',0,'verbose', 0, 'solver', 'sdpnal');

out = optimize(F,P(3,3),options)
double(P)