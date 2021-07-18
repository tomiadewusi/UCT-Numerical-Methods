clear 
clc 
rng(0) 

N = 1000; 
deltat = 1/N; 

S0 = [70;100;90]; 
mu = [0.4; 0.1; 0.12]; 
sigma = [0.4; 0.22; 0.25]; 
Sigma = [1.0  0.3   0.95; 
         0.3  1.0   0.55;
         0.95 0.55  1.00; ]; 
% all(Sigma == Sigma','all')
Spath = S0.*exp(cumsum((mu-0.5.*sigma.^2).*deltat +...
    sigma.*sqrt(deltat).*chol(Sigma,'lower')*randn(3,N),2)); 
% figure() 
% plot(Spath')
% These two should be equivalent they only differ slightly due to machine
% precision 
log_returns = log(Spath(:,2:end)) - log(Spath(:,1:end-1)); 
log_returns1 = log(Spath(:,2:end)./Spath(:,1:end-1)); 
% all(log_returns == log_returns1, 'all') 
% figure() 
% plot(log_returns') 
corrcoef(log_returns')
corrcoef(log_returns1') 