clear 
clc 
rng(0) 

T = 2; 
S0 = 100; 
sigma = 0.45; 
r = 0.1; 
K = 110; 
X = 10; 
 
% Analytical Solution 
d = (log(S0/K) + (r - 0.5*sigma^2)*T)/(sigma*sqrt(T)); 
c = X*exp(-r*T)*normcdf(d); 
delta = (X*exp(-r*T)/(S0*sigma*sqrt(T)))*normpdf(d); 
n =10000; 
Z = randn(n,1); 
ST = S0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*Z);
fx = exp(-r*T)*X*(ST>K);
cMC = mean(fx); 
fxdelta = fx.*Z./(S0*sigma*sqrt(T));
deltaMC = mean(fxdelta); 
confidence_delta = 