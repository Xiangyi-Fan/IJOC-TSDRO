function [epsilon0_sp] =  LS_epsilon_sp_value(data_cons,data)

epsilon_dro_list = [1,2.5,5,10]*20;
K0 = 2;
Npoints = size(data,2);
obj_dro_list = zeros(1,K0);
obj_dro_avg_list = zeros(1,length(epsilon_dro_list));

% K fold preparation
data0 = cell(K0,1);
n = floor(Npoints/K0);
for i = 1:K0
    data0{i,1} = data(:,(i-1)*n+1:i*n);
end


for i = 1:length(epsilon_dro_list)
    epsilon_dro = epsilon_dro_list(i);
    for j = 1:K0
        data_test = data0{j,1};
        len = size(data_test,1);
        data_test = data_test(1:(len-1),:);
        data_train = [];
        for k = 1:K0
            if k ~= j
                data_train = [data_train,data0{k,1}];
            end
        end
        % DRO 
        [x_dro,obj_dro] = LS_PLD_sp_CV(data_cons,data_train,epsilon_dro);
        x_dro = double(x_dro);
        obj_dro_list(j) = obj_exp_LS(x_dro,data_test);       
    end
    obj_dro_avg_list(i) = sum(obj_dro_list)/K0;
end


min_dro = 9999999;
min_dro_index = 0;
for i = 1:length(epsilon_dro_list) 
    if obj_dro_avg_list(i) < min_dro
        min_dro = obj_dro_avg_list(i);
        min_dro_index = i;
    end
end

epsilon0_sp = epsilon_dro_list(min_dro_index);
end