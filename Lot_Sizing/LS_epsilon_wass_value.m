function [epsilon0_wass] =  LS_epsilon_wass_value(data)

epsilon_wass_list = [1,5,10];
K0 = 2;
Npoints = size(data,2);
obj_wass_list = zeros(1,K0);
obj_wass_avg_list = zeros(1,length(epsilon_wass_list));

% K fold preparation
data0 = cell(K0,1);
n = floor(Npoints/K0);
for i = 1:K0
    data0{i,1} = data(:,(i-1)*n+1:i*n);
end


for i = 1:length(epsilon_wass_list)
    epsilon_wass = epsilon_wass_list(i);
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
        
        % wasserstein
        [x_wass,obj_wass] = LS_wass_CV(data_train,epsilon_wass);
        x_wass = double(x_wass);
        obj_wass_list(j) = obj_exp_LS(x_wass,data_test);       
    end
    obj_wass_avg_list(i) = sum(obj_wass_list)/K0;
end


min_wass = 9999999;
min_wass_index = 0;
for i = 1:length(epsilon_wass_list) 
    if obj_wass_avg_list(i) < min_wass
        min_wass = obj_wass_avg_list(i);
        min_wass_index = i;
    end
end

epsilon0_wass = epsilon_wass_list(min_wass_index);
end