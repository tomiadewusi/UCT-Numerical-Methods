clear 
clc 
rng(0) 

n = 5; 
N = 100; 
S0 = 100; 
S100 = 250; 
mu = 0.15; 
sigma = 0.2; 
T = 5; 

t = linspace(0,T,N+1); 
W0 = zeros(n,1); 
W100 = (log(S100/S0) - (mu - 0.5*sigma.^2)*T)/sigma*ones(n,1); 
Z = randn(n,N-1); 

W = [W0 BrownianBridge(t,W0,W100,Z) W100] ; 
ST = S0.*exp((mu - 0.5*sigma.^2).*t + sigma.*W);
figure() 
subplot(1,2,1)
plot(t,W') 
subplot(1,2,2) 
plot(t,ST')