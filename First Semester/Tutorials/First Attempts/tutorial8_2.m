clear 
clc 
rng(0)

T = 2; 
S0 = 100; 
sigma = 0.2; 
r = 0.1; 
K = 100; 

d1 = (log(S0/K) + (r+0.5*sigma^2)*T)/(sigma*sqrt(T)); 
d2 = d1 - sigma*sqrt(T); 
BS = S0*normcdf(d1) - exp(-r*T)*K*normcdf(d2); 

mat = zeros(50,3); 
mat2 = zeros(50,3);
for idx = 1:50 
    n = 1000*idx; 
    Z = randn(n,1);
    ST = S0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T).*Z); 
    fx = exp(-r*T)*max(ST-K ,0); 
    option_value = mean(fx); 
    option_confidence = 3*std(fx)/sqrt(n); 
    lower = BS - option_confidence; 
    upper = BS + option_confidence; 
    mat(idx,:) = [option_value lower upper];
    
    ST_minus = S0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T).*(-Z(1:n/2)));
    fx_minus = exp(-r*T)*max(ST_minus-K ,0);
    fx_dash = exp(-r*T)*max(ST(1:n/2)-K ,0);
    option_value_anth = 0.5*(mean(fx_dash + fx_minus));
    anth_confidence = 1.5*std(fx_dash + fx_minus)/sqrt(n/2); 
    lower_anth = BS - anth_confidence; 
    upper_anth = BS + anth_confidence; 
    mat2(idx,:) = [option_value_anth lower_anth upper_anth]; 
end 

figure()
hold on 
plot([0 50], [BS BS],'k-')
plot(mat(:,1),'b.')
plot(mat(:,[2 3]),'b--')
plot(mat2(:,1),'r.')
plot(mat2(:,[2 3]),'r--')


hold off 