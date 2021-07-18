clear 
clc
rng(0)

n = 5000; 
lambda = 100; 
vector = ncx2rnd(1,lambda, n, 1); 
figure() 
hist(vector) 

p = linspace(0,1,11);
mu = lambda + 1;
sigma = sqrt(2 + 4*lambda);
x = norminv(p,mu,sigma); 

n = sum( vector > x(1:10) & vector < x(2:11) ); 
    
Qd = sum(((n - 500).^2)./500)
lower = chi2inv(0.025, 9)
upper = chi2inv(0.975, 9)

%Since Qd is higher than the upper limit, we reject the null hypothesis 
%if Lambda is 1000 then we fail to reject the null hypothesis 
