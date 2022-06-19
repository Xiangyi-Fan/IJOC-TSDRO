function [A,b,obj,out] = polyhedron_copositive(S,t,dim)
% find (approx) minimum volume ellipsoid covering the
% "set {x: Sx<=t} projected onto first 'dim' dimensions"
% dont input dim if not projecting onto lower dimensions


options = sdpsettings('verbose', 0, 'solver', 'mosek');

[M,N] = size(S);
if nargin == 2
    dim = N;
end
Sc = [-S, t];

A = sdpvar(dim);
b = sdpvar(dim,1);
NN = sdpvar(M);
obj = -logdet(A);

Ac = [A,zeros(dim,N-dim),b];

H = -Sc'*NN*Sc;
H(end,end) = H(end,end) + 1; 

cons = [NN(:) >= 0; [H,Ac';Ac,eye(dim)] >= 0];

out = optimize(cons,obj,options);

A = double(A);b = double(b);
obj = exp(double(obj));
end