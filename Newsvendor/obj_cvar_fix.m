function [obj] =  obj_cvar_fix(x,data,g,Spoints)

N = size(x,1);
delta = 0.1; % risk attitude
s = [10;12];

y1 = zeros(N,Spoints);
y2 = zeros(N,Spoints);
obj_list = zeros(Spoints,1);

for k = 1:Spoints
    xi = data(1:N,k);
    for j = 1:N
        if (x(j) - xi(j)) >= 0
            y1(j,k) = x(j) - xi(j);
            y2(j,k) = 0;
        else
            y1(j,k) = 0;
            y2(j,k) = xi(j) - x(j);
        end
    end
    obj_list(k,1)= g'*y1(:,k) + s'*y2(:,k);      
end

VaR = quantile(obj_list,1-delta);

% calculate CVaR
cvar_list = obj_list(obj_list >= VaR);
obj = mean(cvar_list);

end