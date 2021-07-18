clear
clc
rng(0)

S0 = [35; 65];
sigma = [0.25; 0.2];
Sigma = [1 0.6; 0.6 1];
r = 0.1;
T = 2; 
K = 100; 
L = chol(Sigma,'lower');

n = 50001;
k = ceil(log(n - 1)/log(2)); 

ZS = [Sobol(gendirnums(7,[1 3],k),n);...
    Sobol(gendirnums(11,[1 1 5],k),n)] ; 
ZS(:,1) = []; 

figure()
subplot(1,2,1)
plot(ZS(1,1:1000),ZS(2,1:1000),'k.')
title('Sobol Sequence')
subplot(1,2,2)
plot(norminv(ZS(1,1:1000)),norminv(ZS(2,1:1000)), 'k.')
title('NormInv numbers')

ZS = norminv(ZS); 
mat1 = zeros(50,6); 
for ix = 1:50
    n = ix*1000; 
    
    % Crude Monte Carlo 
    Z_temp = randn(2,n);
    ZN = L*Z_temp;
    
    ST = zeros(2,n);
    ST(1,:) = S0(1)*exp((r-0.5*sigma(1)^2)*T + sigma(1)*sqrt(T).*ZN(1,:));
    ST(2,:) = S0(2)*exp((r-0.5*sigma(2)^2)*T + sigma(2)*sqrt(T).*ZN(2,:));

    sum_ST = sum(ST); 
    
    fx = exp(-r*T)*max(sum_ST - K, 0);
    mat1(ix,1) = mean(fx); 
    mat1(ix,2) = 3*std(fx)/sqrt(n); 
    
    % Pathwise estimate
    dfx = exp(-r*T).*(sum_ST > K).*ST(1,:)./S0(1); 
    mat1(ix,3) = mean(dfx); 
    mat1(ix,4) = 3*std(dfx)/sqrt(n);

    % Sobol Sequence
    Z_temp = ZS(:,1:n);
    ZN = L*Z_temp; 
    
    ST = zeros(2,n);
    ST(1,:) = S0(1)*exp((r-0.5*sigma(1)^2)*T + sigma(1)*sqrt(T).*ZN(1,:));
    ST(2,:) = S0(2)*exp((r-0.5*sigma(2)^2)*T + sigma(2)*sqrt(T).*ZN(2,:));
    
    sum_ST = sum(ST); 
    mat1(ix,5) = mean(exp(-r*T)*max(sum_ST - K, 0));
    
    % Pathwise estimate 
    mat1(ix,6) = mean(exp(-r*T).*(sum_ST > K).*ST(1,:)./S0(1)); 
    
    
end 


best_estimate = mat1(50,5); 
lower = best_estimate - mat1(:,2); 
upper = best_estimate + mat1(:,2); 

best_estimate1 = mat1(50,6); 
lower1 = best_estimate1 - mat1(:,4); 
upper1 = best_estimate1 + mat1(:,4); 


figure()
hold on 
plot( [0 50] , [best_estimate best_estimate], 'k-')
plot( [lower upper], 'b--')
plot(mat1(:,1) ,'b.')
plot(mat1(:,5) ,'r.')
% ylim([21 22.4])
hold off 

figure()
hold on 
plot( [lower1 upper1], 'b--')
plot(mat1(:,3) ,'b.')
plot(mat1(:,6) ,'r.')
% ylim([0.79 0.825])
hold off 
