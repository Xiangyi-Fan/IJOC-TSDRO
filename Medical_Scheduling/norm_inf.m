function [value] =  norm_inf(M0)

M = abs(M0);
value = max(M(:));

end