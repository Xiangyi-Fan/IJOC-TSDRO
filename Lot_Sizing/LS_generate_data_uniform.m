function [data] = LS_generate_data_uniform(Npoints)

N       = 5; % number of locations
v_l = 40; v_u = 50; % v_l <= v <= v_u
u_l = 20 ; u_u = 40; % u_l <= u <= u_u 


data1 = unifrnd(v_l,v_u,[N*N,Npoints]);
v = data1;


data2 = unifrnd(u_l,u_u,[N,Npoints]);
u = data2;


data = [v;u;ones(1,Npoints)];

end