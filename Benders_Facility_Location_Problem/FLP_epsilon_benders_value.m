function [epsilon0] =  FLP_epsilon_benders_value(data,f,C)

epsilon_list = [0.5,5,20];
K0 = 2;
Npoints = size(data,2);
obj_list = zeros(1,K0);
obj_avg_list = zeros(1,length(epsilon_list));

% K fold preparation
data0 = cell(K0,1);
n = floor(Npoints/K0);
for i = 1:K0
    data0{i,1} = data(:,(i-1)*n+1:i*n);
end


for i = 1:length(epsilon_list)
    epsilon_bender = epsilon_list(i);
    for j = 1:K0
        data_test = data0{j,1};
        data_train = [];
        for k = 1:K0
            if k ~= j
                data_train = [data_train,data0{k,1}];
            end
        end
        % copositive 
        [x_benders] = FLP_PLD_benders_CV(data_train,epsilon_bender,f,C);
        obj_list(j) = obj_exp_FLP(data_test,x_bendersf,C);       
    end
    obj_avg_list(i) = sum(obj_list)/K0;
end


min_benders = 9999999;
min_benders_index = 0;
for i = 1:length(epsilon_list) 
    if obj_avg_list(i) < min_benders
        min_benders = obj_avg_list(i);
        min_benders_index = i;
    end
end

epsilon0 = epsilon_list(min_benders_index);
end