function [data] = generate_data(Npoints)

seed = 1;
rng(seed);

N       = 5; % number of stocks

xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u
s_l = 0; s_u = 50; % s_l <= s <= s_u 


mu_list = [0.8;0.9;1;1.1;1.2];
sigma_list = [1;1;1;1;1]*1;

data1 = zeros(N,Npoints);
for i = 1:N
    mu = mu_list(i);
    sigma = sigma_list(i);
    pd2 = makedist('Lognormal','mu',mu,'sigma',sigma);
    t1 = truncate(pd2,xi_l,xi_u);
    data1(i,:) = random(t1,1,Npoints);
end
xi = data1;
    

mu_list = [3.2;3.2;3.3;3.3;3.4];
sigma_list = [1;1;1;1;1]*2;

data2 = zeros(N,Npoints);
for i = 1:N
    mu = mu_list(i);
    sigma = sigma_list(i);
    pd4= makedist('Lognormal','mu',mu,'sigma',sigma);
    t2 = truncate(pd4,s_l,s_u);
    data2(i,:) = random(t2,1,Npoints);
end
s = data2;

data = [xi;s;ones(1,Npoints)];

end