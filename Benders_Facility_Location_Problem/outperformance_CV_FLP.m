function [avg_exp_benders,avg_exp_saa,benders10,benders90,saa10,saa90,avg_T_benders,avg_T_saa,avg_T_max_sub, avg_num_iter] = outperformance_CV_FLP(Npoints,K)

seed = 1;
rng(seed);

trial = 100;
Spoints = 300;

I = 5;
J = 5;
f_list = cell(trial,1);
C_list = cell(trial,1);
for i = 1:trial
    f_list{i} = rand(I,1)*1000+4000;
    C_list{i} = rand(I,J)*90+10;
end

exp_benders_list = zeros(trial,1);
exp_saa_list = zeros(trial,1);

T_benders = zeros(trial,1);
T_saa = zeros(trial,1);

T_sub_benders = zeros(trial,1);
num_iter_benders = zeros(trial,1);

data_out = cell(trial,1);
data_in = cell(trial,1);

for i = 1:trial
    data_out{i,1} = FLP_generate_data(Spoints);
end 

for i = 1:trial
    data_in{i,1} = FLP_generate_data(Npoints);
end

data_cons = FLP_generate_data(K);


disp("finish generating data")



for i = 1:trial
    % Cross Validation
    if Npoints <= 20
        epsilon_benders = 1;
    else
        epsilon_benders = 1;
    end
    
    % copositive programming 
    [benders_optimal] = Benders(data_in{i,1},data_in{i,1},epsilon_benders,K,f_list{i},C_list{i});
    if benders_optimal.iter >= 30
        continue
    end
    x_benders = benders_optimal.q;
    T_benders(i) = benders_optimal.runtime;
    T_sub_benders(i) = benders_optimal.runtime_sub;
    num_iter_benders(i) = benders_optimal.iter;
    
    exp_benders_list(i) = obj_exp_FLP(data_out{i,1},x_benders,f_list{i},C_list{i});
    disp("finish benders")
    
    
    % SAA
    [x_saa,obj_saa,runtime_saa] = FLP_SAA(data_in{i,1},f_list{i},C_list{i});
    T_saa(i) = runtime_saa;
    
    exp_saa_list(i) = obj_exp_FLP(data_out{i,1},x_saa,f_list{i},C_list{i});
    disp("finish SAA")


    X = ['tryyyyyyyyyyyyy ',num2str(i),'th trial for N =', num2str(Npoints),' lalalalalalalala'];
    disp(X)
    
end

exp_benders_list = nonzeros(exp_benders_list);
exp_saa_list = nonzeros(exp_saa_list);
T_benders = nonzeros(T_benders);
T_saa = nonzeros(T_saa);
T_sub_benders = nonzeros(T_sub_benders);
num_iter_benders = nonzeros(num_iter_benders);
avg_exp_benders = mean(exp_benders_list);
benders10 = quantile(exp_benders_list,0.1);
benders90 = quantile(exp_benders_list,0.9);
avg_exp_saa = mean(exp_saa_list);
saa10 = quantile(exp_saa_list,0.1);
saa90 = quantile(exp_saa_list,0.9);
avg_T_benders = mean(T_benders);
avg_T_saa = mean(T_saa);
avg_T_max_sub = mean(T_sub_benders);
avg_num_iter = mean(num_iter_benders);
avg_exp_benders
avg_exp_saa
avg_T_benders
avg_T_saa 
end