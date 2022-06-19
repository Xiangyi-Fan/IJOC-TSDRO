function [data] =  PLD_data(K)

% rng(1);

N       = 3; % number of stocks

xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u
s_l = 20 ; s_u = 40; % s_l <= s <= s_u 

pd1 = makedist('Uniform','Lower',xi_l,'Upper',xi_u);
data1 = random(pd1,N,K);
xi = data1;

pd2 = makedist('Uniform','Lower',s_l,'Upper',s_u);
data2 = random(pd2,N,K);
s = data2;

data = [xi;s;ones(1,K)];

end