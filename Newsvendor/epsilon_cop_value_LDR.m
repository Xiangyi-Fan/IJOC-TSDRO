function [epsilon0_cop] =  epsilon_cop_value_LDR(data)

epsilon_cop_list = [0.5,1,5,10,15,20,30,40,50];
K = 5;
g = [1;3;5];
Npoints = size(data,2);
obj_cop_list = zeros(1,K);
obj_cop_avg_list = zeros(1,length(epsilon_cop_list));

% K fold preparation
data0 = cell(K,1);
n = floor(Npoints/K);
for i = 1:K
    data0{i,1} = data(:,(i-1)*n+1:i*n);
end


for i = 1:length(epsilon_cop_list)
    epsilon_cop = epsilon_cop_list(i);
    for j = 1:K
        data_test = data0{j,1};
        len = size(data_test,1);
        data_test = data_test(1:(len-1),:);
        data_train = [];
        for k = 1:K
            if k ~= j
                data_train = [data_train,data0{k,1}];
            end
        end
        % copositive 
        [x_cop,obj_cop] = Newsvendor_LDR_cop_CV(data_train,epsilon_cop);
        x_cop = double(x_cop);
        obj_cop_list(j) = obj_cvar(x_cop,data_test,g,size(data_test,2));       
    end
    obj_cop_avg_list(i) = sum(obj_cop_list)/K;
end


min_cop = 9999999;
min_cop_index = 0;
for i = 1:length(epsilon_cop_list) 
    if obj_cop_avg_list(i) < min_cop
        min_cop = obj_cop_avg_list(i);
        min_cop_index = i;
    end
end

epsilon0_cop = epsilon_cop_list(min_cop_index);
end