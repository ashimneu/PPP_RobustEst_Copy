function [g1] = lts_g(A,num)
m = size(A,1);
n = size(A,2);

p = m+n;
g1 = floor(p/2) + floor((n+1)/2);
g2 = floor(num/2) + floor((n+1)/2);
end

