seed = 1;
rng(seed);

trial = 100;
N = 2; % number of stocks
Spoints = 10000;
Npoints = 20;
xi_l = 0; xi_u = 10; % xi_l <= xi <= xi_u
g = [2;4];
s = [10;12];
   
obj_dro_list = zeros(trial,1);
obj_saa_list = zeros(trial,1);
cvar_dro_list = zeros(trial,1);
cvar_saa_list = zeros(trial,1);
gap_saa = zeros(trial,1);

data_out = cell(trial,1);
data_in = cell(trial,1);

for i = 1:trial
    
    data_out{i,1} = generate_data_fix(Spoints);
    
    data_in{i,1} = generate_data_fix(Npoints);
end

disp("finish generating data")

for i = 1:trial
    
    % copositive programming    
    [x_dro,obj_dro] = Newsvendor_PLD_cop_fix(data_in{i,1});
    x_dro = double(x_dro);
    obj_dro = double(obj_dro);
    
    cvar_dro_list(i) = obj_cvar_fix(x_dro,data_out{i,1},g,Spoints);
    disp("finish DRO")
    
    
    % SAA
    [x_saa,obj_saa] = SAA_fix(data_in{i,1});
    x_saa = double(x_saa);
    obj_saa = double(obj_saa);
    
    cvar_saa_list(i) = obj_cvar_fix(x_saa,data_out{i,1},g,Spoints);
    disp("finish SAA")
    
    gap_saa(i) = (cvar_saa_list(i) - cvar_dro_list(i))/cvar_saa_list(i);
    disp("tryyyyyyyyyyyyy next trial lalalalalalalala")
end

avg_cvar_dro = mean(cvar_dro_list);
avg_cvar_saa = mean(cvar_saa_list);
avg_gap_saa = mean(gap_saa);
avg_cvar_dro
avg_cvar_saa
avg_gap_saa