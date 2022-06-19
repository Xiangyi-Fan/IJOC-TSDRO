function [P_k] =  FLP_PLD_partitions(data,k)

J       = 5; % number of stocks
K = size(data,2);

d_l =50; d_u = 300; % v_l <= v <= v_u
ee1 = ones(J,1);
P = [eye(J), -1*d_l*ee1;
     -1*eye(J), d_u*ee1];

data1 = data(1:(size(data,1) - 1),:);
for i = 1:K
    if i ~= k
        P = [P; [-2*(data1(:,i)-data1(:,k))',data(:,i)'*data(:,i)-data(:,k)'*data(:,k)]];
    end
end

P_k = P;

end
