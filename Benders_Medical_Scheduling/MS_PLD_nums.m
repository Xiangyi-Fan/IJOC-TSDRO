function [I_k] =  MS_PLD_nums(P_k, data)
I_k = 0;
Npoints = size(data,2); % number of data points

for i = 1:Npoints
    if P_k*data(:,i) >= 0
        I_k = I_k + 1;
    end
end

end