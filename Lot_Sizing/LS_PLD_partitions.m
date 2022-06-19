function [P_k] =  LS_PLD_partitions(data,k)

N       = 5; % number of stocks
K = size(data,2);

v_l = 40; v_u = 50; % v_l <= v <= v_u
u_l = 20 ; u_u = 40; % u_l <= u <= u_u 
ee1 = ones(N*N,1);
ee2 = ones(N,1);
P = [eye(N*N), zeros(N*N,N), -1*v_l*ee1;
     -1*eye(N*N), zeros(N*N,N), v_u*ee1;
     zeros(N,N*N), eye(N), -1*u_l*ee2;
     zeros(N,N*N), -1*eye(N), u_u*ee2];

data1 = data(1:(size(data,1) - 1),:);
for i = 1:K
    if i ~= k
        P = [P; [-2*(data1(:,i)-data1(:,k))',data(:,i)'*data(:,i)-data(:,k)'*data(:,k)]];
    end
end

P_k = P;

end
