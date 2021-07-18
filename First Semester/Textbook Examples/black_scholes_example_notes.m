clear 
clc 
rng(0)

T = 2; 
r = 0.085; 
sigma = 0.3;
K = 100; 
S0 = 100; 

% Closed Form Solution
actual = S0*normcdf( (log(S0/K)+(r + 0.5*sigma.^2)*T )/(sigma*sqrt(T)) )...
    -exp(-r*T)*K* normcdf(log(S0/K)+(r - 0.5*sigma.^2)*T/(sigma*sqrt(T)) ); 

t = 1:50; 
mat = zeros(length(t),4); 
for idx=t
    N = idx*1000;
    Z = randn(N,1);
    f_Z = exp(-r*T).*(max(S0*exp((r-0.5*sigma.^2)*T + sigma*sqrt(T).*Z)-K,0));
    value = sum(f_Z)/N; 
    confidence = 3*(sqrt(sum((f_Z -value).^2)/(N-1)))/sqrt(N);
    upper = actual + confidence; 
    lower = actual - confidence; 
    mat(idx,:) = [value actual upper lower]; 
end 

figure()
hold on
plot(t, mat(:,2:end) )
plot(t,mat(:,1),'k.')

