function [P_k] =  PLD_partitions_fix(data,k)

N       = 2; % number of stocks
K = size(data,2);

xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u
ee = ones(N,1);
P = [eye(N), -1*xi_l*ee;
     -1*eye(N), xi_u*ee];

data1 = data(1:(size(data,1)-1),:);
for i = 1:K
    if i ~= k
        P = [P; [-2*(data1(:,i)-data1(:,k))',data(:,i)'*data(:,i)-data(:,k)'*data(:,k)]];
    end
end

P_k = P;

end