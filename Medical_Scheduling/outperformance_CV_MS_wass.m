function [avg_cvar_cop,avg_cvar_sp,avg_cvar_saa,avg_cvar_wass,cop10,cop90,sp10,sp90,saa10,saa90,wass10,wass90,avg_T_cop,avg_T_sp,avg_T_saa,avg_T_wass] = outperformance_CV_MS_wass(Npoints,K,epsilon_cop_fix,epsilon_sp_fix)

seed = 1;
rng(seed);

trial = 50;
N = 8; % number of patients
Spoints = 10000;
   

cvar_cop_list = zeros(trial,1);
cvar_sp_list = zeros(trial,1);
cvar_saa_list = zeros(trial,1);
cvar_wass_list = zeros(trial,1);
gap_saa = zeros(trial,1);
gap_sp = zeros(trial,1);
gap_wass = zeros(trial,1);

T_cop = zeros(trial,1);
T_sp = zeros(trial,1);
T_saa = zeros(trial,1);
T_wass = zeros(trial,1);

data_out = cell(trial,1);
data_in = cell(trial,1);

for i = 1:trial
    data_out{i,1} = MS_generate_data(Spoints);
end 

for i = 1:trial
    data_in{i,1} = MS_generate_data(Npoints);
end

data_cons =  MS_generate_data(K);

seed = 1;
rng(seed);
for i = 1:trial
    data_out{i,1} = MS_generate_data(Spoints);
end 

disp("finish generating data")

cop_cont = true;
sp_cont = true;
if Npoints > 80
    cop_cont = false;
end


i = 0;
while i < trial
    i = i+1;
    
    % Cross Validation for epsilon

    epsilon_cop = epsilon_cop_fix;
    epsilon_sp = epsilon_sp_fix;
    
    % copositive programming 
    if cop_cont
        [x_cop,obj_cop,runtime_cop] = MS_PLD_cop(data_cons,data_in{i,1},epsilon_cop,K);
        if runtime_cop > 30*60
            cop_cont = false;
            T_cop(i) = NaN;
            cvar_cop_list(i) = NaN;
        else  
            T_cop(i) = runtime_cop;
            x_cop = double(x_cop);
            cvar_cop_list(i) = obj_cvar_MS(x_cop,data_out{i,1});
        end
    else
        T_cop(i) = NaN;
        cvar_cop_list(i) = NaN;       
    end
    disp("finish copositive")
    
    
    % S procedure programming    
    if sp_cont
        [x_sp,obj_sp,runtime_sp] = MS_PLD_sp(data_cons,data_in{i,1},epsilon_sp,K);
        if runtime_sp > 30*60
            sp_cont = false;
            T_sp(i) = NaN;
            cvar_sp_list(i) = NaN;
        else  
            T_sp(i) = runtime_sp;
            x_sp = double(x_sp);
            cvar_sp_list(i) = obj_cvar_MS(x_sp,data_out{i,1});
        end
    else
        T_sp(i) = NaN;
        cvar_sp_list(i) = NaN;       
    end
    disp("finish S lemma")
    
    
    % SAA
    
    [x_saa,obj_saa,runtime_saa] = MS_cvar_SAA(data_in{i,1});
    T_saa(i) = runtime_saa;
    x_saa = double(x_saa);
    obj_saa = double(obj_saa);
    
    cvar_saa_list(i) = obj_cvar_MS(x_saa,data_out{i,1});
    disp("finish SAA")
   

    gap_saa(i) = (cvar_saa_list(i) - cvar_cop_list(i))/cvar_saa_list(i);
    gap_sp(i) = (cvar_sp_list(i) - cvar_cop_list(i))/cvar_sp_list(i);
    X = ['tryyyyyyyyyyyyy',num2str(i),'th trial for N =', num2str(Npoints),'lalalalalalalala'];
    disp(X)
    
end


avg_cvar_cop = mean(cvar_cop_list);
cop10 = quantile(cvar_cop_list,0.1);
cop90 = quantile(cvar_cop_list,0.9);
avg_cvar_sp = mean(cvar_sp_list);
sp10 = quantile(cvar_sp_list,0.1);
sp90 = quantile(cvar_sp_list,0.9);
avg_cvar_saa = mean(cvar_saa_list);
saa10 = quantile(cvar_saa_list,0.1);
saa90 = quantile(cvar_saa_list,0.9);
avg_cvar_wass = mean(cvar_wass_list);
wass10 = quantile(cvar_wass_list,0.1);
wass90 = quantile(cvar_wass_list,0.9);
avg_T_cop = mean(T_cop);
avg_T_sp = mean(T_sp);
avg_T_saa = mean(T_saa);
avg_T_wass = mean(T_wass);
avg_gap_saa = mean(gap_saa);
avg_gap_sp = mean(gap_sp);
avg_gap_wass = mean(gap_wass);
avg_cvar_cop
avg_cvar_sp
avg_cvar_saa
avg_cvar_wass
avg_gap_saa
avg_gap_sp
avg_gap_wass
avg_T_cop
avg_T_sp
avg_T_saa
avg_T_wass
end