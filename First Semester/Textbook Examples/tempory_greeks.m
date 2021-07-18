clear
clc
rng(0)

S0 = 100;
r = 0.085;
sigma = 0.3;
T = 2;
K = 100;
delta_S0 = 5;
S0_1 = S0 + delta_S0;
S0_2 = S0;
 
d1 = (log(S0/K) + T*(r + 0.5*sigma.^2) )/(sigma*sqrt(T));  
d2 = d1 - sigma*sqrt(T); 
act_call_delta = normcdf(d1,0,1)
