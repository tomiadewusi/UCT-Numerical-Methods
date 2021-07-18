clear 
clc 
rng(0)

s01 = 100; 
s02 = 90; 
mu1 = 0.2; 
mu2 = 0.3;
sigma1 = 0.25; 
sigma2 = 0.45; 
p = 0.9; 
t = 0.25; 

Z = randn(2,1000000);
smatrix = [1 p; p 1]; 
L = chol(smatrix,'lower'); 
%To note that the convention in the notes is to use a lower triangular
%matrix 

X = L*Z; 

S = zeros(size(X)); 
S(1,:) = s01.*exp((mu1 - 0.5.*sigma1.^2).*t + sigma1.*sqrt(t).*X(1,:));
S(2,:) = s02.*exp((mu2 - 0.5.*sigma2.^2).*t + sigma2.*sqrt(t).*X(2,:));

figure(1) 
subplot(1,2,1) 
hist(S(1,:), 50) 
subplot(1,2,2)
hist(S(2,:), 50) 

figure(2) 
plot(S(1,:),S(2,:),'.')

figure(3)
hist3(S',[50 50])

%Sample mean and variance 
m1 = mean(S(1,:)); 
m2 = mean(S(2,:)); 
v1 = var(S(1,:)); 
v2 = var(S(2,:));

%theorectical mean and variance 
mt1 = s01*exp(mu1*t); 
mt2 = s02*exp(mu2*t); 
vt1 = (mt1^2)*(exp(sigma1^2*t)-1); 
vt2 = (mt2^2)*(exp(sigma2^2*t)-1); 



