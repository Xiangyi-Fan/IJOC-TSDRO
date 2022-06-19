function [data] = LS_generate_data(Npoints)

N       = 5; % number of locations
v_l = 40; v_u = 50; % v_l <= v <= v_u
u_l = 20 ; u_u = 40; % u_l <= u <= u_u 


mu_list = ones(N*N,1)*3; 
sigma_list = ones(N*N,1)*0.1;%(first 0.5, second 0.1,third 0.2)

data1 = zeros(N*N,Npoints);
for i = 1:N*N
    mu = mu_list(i);
    sigma = sigma_list(i);
    pd2 = makedist('Lognormal','mu',mu,'sigma',sigma);
    t1 = truncate(pd2,v_l,v_u);
    data1(i,:) = random(t1,1,Npoints);
end
v = data1;

mu_list = [3;3;3;3.5;3.5;3.5];%(first [3;3;3;3;3;3], second [3;3;3;4;4;4])
sigma_list = ones(N,1)*0.2;%(first 0.2, second 0.2,third 0.1)

data2 = zeros(N,Npoints);
for i = 1:N
    mu = mu_list(i);
    sigma = sigma_list(i);
    pd4= makedist('Lognormal','mu',mu,'sigma',sigma);
    t2 = truncate(pd4,u_l,u_u);
    data2(i,:) = random(t2,1,Npoints);
end
u = data2;

data = [v;u;ones(1,Npoints)];

end