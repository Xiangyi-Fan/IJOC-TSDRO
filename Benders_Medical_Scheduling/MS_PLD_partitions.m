function [P_k] =  MS_PLD_partitions(data,k,K)

N       = 8; % number of patients

l_l = 20; l_u = 100; % l_l <= l <= l_u
pi_l = 1 ; pi_u = 10; % pi_l <= pi <= pi_u 

ee = ones(N,1);
P = [eye(N), zeros(N), -1*l_l*ee;
     -1*eye(N), zeros(N), l_u*ee;
     zeros(N), eye(N), -1*pi_l*ee;
     zeros(N), -1*eye(N), pi_u*ee];

data1 = data(1:(size(data,1) - 1),:);
for i = 1:K
    if i ~= k
        P = [P; [-2*(data1(:,i)-data1(:,k))',data(:,i)'*data(:,i)-data(:,k)'*data(:,k)]];
    end
end

P_k = P;

end
