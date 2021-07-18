clear 
clc
rng(0)

Z = randn(5,99); 
T = 5; 
t = 0:0.05:T;
S0 = 100; 
S100 = 250; 
mu = 0.15; 
sigma = 0.2; 
W0 = zeros(5,1); 
W100 = ((log(S100/S0) - (mu - 0.5*sigma^2)*T)/sigma).*ones(5,1); 
W = [W0 BrownianBridgeV(t,W0,W100,Z ) W100];  
S = S0*exp((mu - 0.5*sigma^2)*t + sigma.*W);

figure()
subplot(1,2,1)
plot(t,W.')
subplot(1,2,2)
plot(t,S.')