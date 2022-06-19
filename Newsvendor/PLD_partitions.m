function [P_k] =  PLD_partitions(data,k,K)

N       = 5; % number of stocks
% K = size(data,2);

xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u
s_l = 0; s_u = 50; % s_l <= s <= s_u 
ee = ones(N,1);
P = [eye(N), zeros(N), -1*xi_l*ee;
     -1*eye(N), zeros(N), xi_u*ee;
     zeros(N), eye(N), -1*s_l*ee;
     zeros(N), -1*eye(N), s_u*ee];

% data = generate_data(K);
data1 = data(1:(size(data,1) - 1),:);

for i = 1:K
    if i ~= k
        P = [P; [-2*(data1(:,i)-data1(:,k))',data(:,i)'*data(:,i)-data(:,k)'*data(:,k)]];
    end
end

P_k = P;

end
