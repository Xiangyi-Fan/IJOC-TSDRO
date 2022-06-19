function [avg_cop,avg_sp,avg_saa,avg_wass,cop0,cop100,sp0,sp100,saa0,saa100,wass0,wass100,feasibility_num,avg_T_cop,avg_T_sp,avg_T_saa,avg_T_wass] = outperformance_LS_general_thm(Npoints)

seed = 1;
rng(seed);

trial = 10;
trial_new = trial*1;
N = 5; % number of stocks
Spoints = 20000;
K = Npoints;
rho1 = 0.1;

dro_cop_list = zeros(trial_new,1);
dro_sp_list = zeros(trial_new,1);
saa_list = zeros(trial_new,1);
wass_list = zeros(trial_new,1);
gap_saa = zeros(trial_new,1);
gap_wass = zeros(trial_new,1);

T_dro_cop = zeros(trial_new,1);
T_dro_sp = zeros(trial_new,1);
T_saa = zeros(trial_new,1);
T_wass = zeros(trial_new,1);

data_out = cell(trial_new,1);
data_in = cell(trial_new,1);
data0_in = cell(trial_new,1);

epsilon = cell(K,1);

for i = 1:trial_new
    data_out{i,1} = LS_generate_data(Spoints);
end 

for i = 1:trial_new
    data_in{i,1} = LS_generate_data(Npoints);
end


seed = 1;
rng(seed);
for i = 1:trial_new
    data_out{i,1} = LS_generate_data(Spoints);
end 
for i = 1:trial_new
    data0_in{i,1} = LS_generate_data(20);
end

disp("finish generating data")


cop_cont = true;
sp_cont = true;
wass_cont = false;
if Npoints > 20
    wass_cont = false;
end
if Npoints > 80
    cop_cont = false;
end

feasibility_num = 0;
for i = 1:trial_new
    
    cons_data = data_in{i,:};
    for k = 1:K
        epsilon_k = LS_PLD_epsilon_thm(rho1, data_in{i,:}, cons_data, k);
        epsilon{k,1} = epsilon_k;
    end  
    
    % coposiitive 
    if cop_cont
        [x_cop,obj_cop,runtime_cop] = LS_PLD_cop_general_thm(data_in{i,1},epsilon);   
        if runtime_cop > 15*60
            cop_cont = false;
            T_dro_cop(i) = NaN;
            dro_cop_list(i) = NaN;
        else       
            T_dro_cop(i) = runtime_cop;
            x_cop = double(x_cop);
            dro_cop_list(i) = obj_exp_LS(x_cop,data_out{i,1});
        end       
    else
        T_dro_cop(i) = NaN;
        dro_cop_list(i) = NaN;
    end
    disp("finish copositive")

    % Slemma 
    if sp_cont
        [x_sp,obj_sp,runtime_sp] = LS_PLD_sp_general_thm(data_in{i,1},epsilon);   
        if runtime_sp > 30*60
            sp_cont = false;
            T_dro_sp(i) = NaN;
            dro_sp_list(i) = NaN;
        else       
            T_dro_sp(i) = runtime_sp;
            x_sp = double(x_sp);
            dro_sp_list(i) = obj_exp_LS(x_sp,data_out{i,1});
        end       
    else
        T_dro_sp(i) = NaN;
        dro_sp_list(i) = NaN;
    end
    disp("finish Slemma")

    % SAA
    
    [x_saa,obj_saa,runtime_saa] = LS_SAA(data_in{i,1});
    T_saa(i) = runtime_saa;
    x_saa = double(x_saa);
    obj_saa = double(obj_saa);
    
    [obj_saa_out,feasibility] = obj_exp_LS_SAA(x_saa,data_out{i,1});
    if feasibility == 0
        feasibility_num = feasibility_num + 1;
        saa_list(feasibility_num) = obj_saa_out;
    end
    disp("finish SAA")

    gap_saa(i) = (saa_list(i) - dro_cop_list(i))/saa_list(i);  
    X = ['tryyyyyyyyyyyyy',num2str(i),'th trial for N =', num2str(Npoints),'lalalalalalalala'];
    disp(X)
end


avg_cop = mean(dro_cop_list);
cop0 = quantile(dro_cop_list,0);
cop100 = quantile(dro_cop_list,1);
avg_sp = mean(dro_sp_list);
sp0 = quantile(dro_sp_list,0);
sp100 = quantile(dro_sp_list,1);
saa_list_nonzero = nonzeros(saa_list);
avg_saa = mean(saa_list_nonzero);
saa0 = quantile(saa_list_nonzero,0);
saa100 = quantile(saa_list_nonzero,1);
avg_wass = mean(wass_list);
wass0 = quantile(wass_list,0);
wass100 = quantile(wass_list,1);
avg_gap_saa = mean(gap_saa);
avg_gap_wass = mean(gap_wass);
avg_T_cop = mean(T_dro_cop);
avg_T_sp = mean(T_dro_sp);
avg_T_saa = mean(T_saa);
avg_T_wass = mean(T_wass);
avg_cop
avg_sp
avg_saa
avg_wass
avg_gap_saa
avg_gap_wass
feasibility_num
avg_T_cop
avg_T_sp
avg_T_saa
avg_T_wass
end
