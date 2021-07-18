clear 
clc
rng(0) 

% Part a

deltat = 2; 
ti = 0:deltat:20; 
N = length(ti) - 1; 
M = N; 
sigma = 0.25; 
n  = 50000; 

% Observed Market Prices 
B0t = exp( (-7.*ti + (t+1).*log(t+1))./(100) ); 
figure() 
plot(ti,B0t) 