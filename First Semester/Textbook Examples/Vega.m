clear
clc
% rng(0)

% Fixed Inputs
S0 = 100;
K = 80;
r = 0.1;
sigma = 0.1;
T = 2; % This is the maturity time
eta = 1;
Bump = 0.000001;
d1 = (log(S0/K) + (r + 0.5*sigma^2).*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
VegaExact = S0*normpdf(d1)*sqrt(T);

ST = @(Z) S0*exp((r-0.5*sigma.^2)*T + sigma*sqrt(T).*Z)  ;

g = @(Z,sigma) exp(-r*T)*max(S0*exp((r-0.5*sigma.^2)*T...
    + sigma*sqrt(T).*Z) - K,0) ;
fpath = @(Z) (exp(-r*T).*(ST(Z) > K)...
        .*ST(Z).*(-sigma*T + sqrt(T).*Z)); 
mat = zeros(200,2);
mat1 = zeros(200,2); 

for ix = 1:200
    n = 1000*ix; 
    Z = randn(n,1);
    % I think that this is fine
    path = fpath(Z); 
    mat(ix,:) = [mean(path) 3*std(path)/sqrt(n)];
    hx = g(Z,sigma + Bump) - g(Z,sigma-Bump)/(2*Bump);  
    mat1(ix,:) = [mean(g(Z,sigma + Bump) - g(Z,sigma-Bump))/(2*Bump)...
        3*std(hx)/sqrt(n)];
    
end
lower = VegaExact - mat(:,2); 
upper = VegaExact + mat(:,2); 
lower1 = VegaExact - mat1(:,2); 
upper1 = VegaExact + mat1(:,2); 

figure() 
subplot(1,2,1) 
hold on 
plot([0 200], [VegaExact VegaExact],'k-')
plot(mat(:,1),'b.') 
plot(lower, 'b--') 
plot(upper, 'b--') 
title('Pathwise Estimate') 
hold off 

subplot(1,2,2) 
hold on 
plot([0 200], [VegaExact VegaExact],'k-')
plot(mat1(:,1),'b.') 
% plot(lower1, 'b--') 
% plot(upper1, 'b--') 
title('Central Difference estimate') 
hold off 


