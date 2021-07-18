clear 
clc 
rng(0) 

S0 = [100;90]; 
mu = [0.2;0.3]; 
sigma = [0.24;0.45]; 
rho = 0.9; 
Sigma = [1 rho; rho 1]; 
t = 0.25; 
Z = randn(2,1000000) ;
X = chol(Sigma,'lower')*Z; 
ST = S0.*exp((mu - 0.5*sigma.^2).*t - sigma*sqrt(t).*X); 
figure() 
subplot(1,2,1) 
hist(ST(1,:),50) 
subplot(1,2,2)
hist(ST(2,:),50) 
figure()
scatter(ST(1,:), ST(2,:),'k.' )
figure() 
hist3(ST', [50 50])



