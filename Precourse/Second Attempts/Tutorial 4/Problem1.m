clear 
clc

f = @(x) x.^3 + 5.*x.^2 + x - 3; 
a = 0; 
b = 2; 
tol = 1e-6;

format long 
value = bisect(f,a,b,tol); 