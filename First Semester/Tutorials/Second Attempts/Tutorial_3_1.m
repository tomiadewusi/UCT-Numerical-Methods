clear 
clc 
rng(0) 
n = 5000; 
lambda = 100; 
RV = ncx2rnd(1,lambda,n,1); 

% figure() 
% hist(RV) 
d = 9; 
p = linspace(0,1,d+2); 
mu = lambda + 1; 
sigma = sqrt(2 + 4*lambda); 
x = norminv(p,mu,sigma); 
Enumpoints = n*p(2); 
Anumpoints = sum((RV < x(2:end)) & ((x(1:end-1) < RV))); 
Qd = sum((Anumpoints - Enumpoints).^2)/Enumpoints;  
alpha = 0.05; 
values = [alpha/2 1-alpha/2]; 
testvalues = chi2inv(values, d)
Qd
  
% We reject the hypothesis 
rng(0) 
n = 5000; 
lambda = 1000; 
RV = ncx2rnd(1,lambda,n,1); 

d = 9; 
p = linspace(0,1,d+2); 
mu = lambda + 1; 
sigma = sqrt(2 + 4*lambda); 
x = norminv(p,mu,sigma); 
Enumpoints = n*p(2); 
Anumpoints = sum((RV < x(2:end)) & ((x(1:end-1) < RV))); 
Qd = sum((Anumpoints - Enumpoints).^2)/Enumpoints;  
alpha = 0.05; 
values = [alpha/2 1-alpha/2]; 
testvalues = chi2inv(values, d)
Qd

% We fail ot reject the hypothesis 