clear 
clc 

mu = 0.3; 
sigma = 0.1; 
x0 = 10; 
delta_t = 0.01; 
t = 0:delta_t:5; 
k = 10; %this is the number of realisations 

num = randn(length(t) - 1,k); 

temp = x0.*exp(cumsum( (mu - 0.5.*sigma.^2).*delta_t + ...
sigma.*sqrt(delta_t).*num));

X = [repmat(x0,1,k) ; temp ]; 

figure(1)
hold on 
plot( repmat(t',1,k), X)
plot(t, x0.*exp(mu.*t), 'k-')
hold off 
