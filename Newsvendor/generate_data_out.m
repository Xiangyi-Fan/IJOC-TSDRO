function [data] =  generate_data_out(Spoints)

N       = 3; % number of stocks

xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u
s_l = 20 ; s_u = 40; % s_l <= s <= s_u 
    
    
    pd1 = makedist('Uniform','Lower',0,'Upper',2);
    mu_list = random(pd1,N,1);
    sigma_list = [1/4;1/4;1/4];
    data1 = zeros(N,Spoints);
    for i = 1:N
        mu = mu_list(i);
        sigma = sigma_list(i);
        pd2 = makedist('Lognormal','mu',mu,'sigma',sigma);
        t1 = truncate(pd2,xi_l,xi_u);
        data1(i,:) = random(t1,1,Spoints);
    end
    xi = data1;

    pd3 = makedist('Uniform','Lower',3.2,'Upper',3.3);
    mu_list = random(pd3,N,1);
    sigma_list = [1/4;1/4;1/4];
    data2 = zeros(N,Spoints);
    for i = 1:N
        mu = mu_list(i);
        sigma = sigma_list(i);
        pd4= makedist('Lognormal','mu',mu,'sigma',sigma);
        t2 = truncate(pd4,s_l,s_u);
        data2(i,:) = random(t2,1,Spoints);
    end
    s = data2;

    data = [xi;s];

end