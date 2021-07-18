clear
clc
rng(0)

T = 2;
S0 = 100;
sigma = 0.4;
r = 0.1;
K = 110;
X = 10;

d = ( log(S0/K) + (r - 0.5*sigma.^2)*T )/(sigma*sqrt(T));
c = X*exp(-r*T)*normcdf(d);
delta_exact = (normpdf(d)*X*exp(-r*T))/(S0*sigma*sqrt(T)); 

mat1 = zeros(50,3);
mat2 = zeros(50,3);
for idx = 1:50
    n = 1000*idx;
    Z = randn(n,1);
    ST = S0*exp((r-0.5*sigma.^2)*T + sigma*sqrt(T).*Z);
    indicator = ST > K;
    payoff = exp(-r*T)*X.*indicator;
    option_value = mean(payoff);
    
    delta  = exp(-r*T)*(Z.*X.*indicator)./(S0*sigma*sqrt(T)); 
    delta_value = mean(delta); 
    conf_delta = 3*std(delta)/sqrt(n); 
    upper_delta = delta_exact + conf_delta; 
    lower_delta = delta_exact - conf_delta; 
  
    confidence = 3*std(payoff)/sqrt(n);
    upper = c + confidence;
    lower = c - confidence;
    mat1(idx,:) = [option_value upper lower];
    mat2(idx,:) = [delta_value upper_delta lower_delta]; 
    
end

figure()
subplot(1,2,1)
hold on 
plot([0 50], [c c], 'b-')
plot(mat1(:,1), 'k.')
plot(mat1(:,[2 3]),'b--')
title('Option value')
%axis tight 
hold off 

subplot(1,2,2)
hold on 
plot([0 50], [delta_exact delta_exact], 'b-')
plot(mat2(:,1), 'k.')
plot(mat2(:,[2 3]),'b--')
title('Delta value using Likelihood estimate') 
%axis tight 
hold off 