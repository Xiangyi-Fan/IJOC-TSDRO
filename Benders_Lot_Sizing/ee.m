function [vec] = ee(J,j)
    e0 = zeros(J,1);
    e0(j) = 1;
    vec = e0;
end