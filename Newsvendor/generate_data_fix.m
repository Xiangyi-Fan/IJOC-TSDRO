function [data] =  generate_data_fix(Npoints)

N       = 2; % number of stocks
xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u


pd1 = makedist('Uniform','Lower',0,'Upper',2);
mu_list = [1;1.5];
sigma_list = [1/4;1/4];
data1 = zeros(N,Npoints);
for i = 1:N
    mu = mu_list(i);
    sigma = sigma_list(i);
    pd2 = makedist('Lognormal','mu',mu,'sigma',sigma);
    t1 = truncate(pd2,xi_l,xi_u);
    data1(i,:) = random(t1,1,Npoints);
end
xi = data1;


data = xi;

end