function [avg_T_benders, avg_T_max_sub, avg_num_iter] = outperformance_LS(Npoints,K)

seed = 1;
rng(seed);

trial = 50;
epsilon0 = 10;

T_benders = zeros(trial,1);
T_sub_benders = zeros(trial,1);
num_iter_benders = zeros(trial,1);
data_in = cell(trial,1);
for i = 1:trial
    data_in{i,1} = LS_generate_data(Npoints);
end


disp("finish generating data")


for i = 1:trial
  
    optimal = Benders(data_in{i,1},data_in{i,1},epsilon0,K);
    T_benders(i) = optimal.runtime;
    T_sub_benders(i) = optimal.runtime_sub;
    num_iter_benders(i) = optimal.iter;  
 
    X = ['Finish ',num2str(i),'th trial for N =', num2str(Npoints)];
    disp(X)    
end

T_benders = nonzeros(T_benders);
T_sub_benders = nonzeros(T_sub_benders);
num_iter_benders = nonzeros(num_iter_benders);
avg_T_benders = mean(T_benders);
avg_T_max_sub = mean(T_sub_benders);
avg_num_iter = mean(num_iter_benders);

avg_T_benders

end