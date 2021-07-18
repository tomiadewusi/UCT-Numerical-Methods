clear 
clc

f = @(x) x^3 + 5*x^2 + x - 3; 

value = bisect(f,0,2,1e-6);