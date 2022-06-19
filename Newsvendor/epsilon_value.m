function [epsilon0_cop, epsilon0_sp] =  epsilon_value(data)

epsilon_cop_list = [0.01,0.1,0.5,1,1.5,2,2.5,3];
epsilon_sp_list = [0.01,0.1,0.5,1,1.5,2,2.5,3];
K = 5;
g = [1;3;5];
Npoints = size(data,2);
obj_sp_list = zeros(1,K);
obj_sp_avg_list = zeros(1,length(epsilon_cop_list));
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
    epsilon_sp = epsilon_sp_list(i);
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
        [x_cop,obj_cop] = Newsvendor_PLD_cop_CV(data_train,epsilon_cop);
        x_cop = double(x_cop);
        obj_cop_list(j) = obj_cvar(x_cop,data_test,g,size(data_test,2));
        
        % S-procedure
        [x_sp,obj_sp] = Newsvendor_PLD_sp_CV(data_train,epsilon_sp);
        x_sp = double(x_sp);
        obj_sp_list(j) = obj_cvar(x_sp,data_test,g,size(data_test,2));        
    end
    obj_cop_avg_list(i) = sum(obj_cop_list)/K;
    obj_sp_avg_list(i) = sum(obj_sp_list)/K;
end


min_cop = 9999999;
min_cop_index = 0;
min_sp = 9999999;
min_sp_index = 0; 
for i = 1:length(epsilon_cop_list) 
    if obj_cop_avg_list(i) < min_cop
        min_cop = obj_cop_avg_list(i);
        min_cop_index = i;
    end
    if obj_sp_avg_list(i) < min_sp
        min_sp = obj_sp_avg_list(i);
        min_sp_index = i;
    end
end

epsilon0_cop = epsilon_cop_list(min_cop_index);
epsilon0_sp = epsilon_sp_list(min_sp_index);
