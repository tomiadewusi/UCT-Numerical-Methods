clear 
clc 
rng(0) 

clear 
clc 
rng(0) 

S0 = 50;
K = 55;
sigma = 0.25;
r = 0.1;
T = 2;
t = 0;
lambda = 0.125; 
delta = 0.2; 

d1 = ((log(S0/K)+(r+0.5*sigma.^2)*(T-t))...
    /(sigma*sqrt(T-t)));
d2 = ((log(S0/K)+(r-0.5*sigma.^2)*(T-t))...
    /(sigma*sqrt(T-t)));
X0 = S0*normcdf(d1) - K*exp(-r*(T-t))*normcdf(d2);

BS = @(S0,K,sigma,r,T,t) ...
    S0.*normcdf(((log(S0./K)+(r+0.5*sigma.^2).*(T-t))...
    ./(sigma.*sqrt(T-t))))...
    - K.*exp(-r.*(T-t)).*normcdf(((log(S0./K)+(r-0.5.*sigma.^2).*(T-t))...
    ./(sigma.*sqrt(T-t))) );

mat1 = zeros(50,2);
mat2 = zeros(50,2);
for ix = 1:50
    n = 1000*ix ; 
    U = rand(n,1) ; 
    tau = -log(U)/lambda; 
    Indicator = tau <= T; 
    mat1(ix,:) = [mean(Indicator), 3*std(Indicator)/sqrt(n)] ; 
    
    tTilde = min(tau,T); 
    Stilde = S0*exp((r-0.5*sigma.^2)*(tTilde)...
        +sigma*sqrt(tTilde).*randn(n,1));
    CVA = X0 - exp(-r.*tTilde).*(max(Stilde-K,0).*(1-Indicator)...
        + delta.*BS(Stilde,K,sigma,r,T,tTilde).*Indicator); 
    mat2(ix,:) = [mean(CVA), 3*std(CVA)/sqrt(n)] ; 
end 

PD = 1 - exp(-lambda*T);
AnaCVA = BS(S0,K,sigma,r,T,0).*(1 - exp(-lambda*T))*(1-delta); 

tauLower = PD - mat1(:,2); 
tauUpper = PD + mat1(:,2); 
CVALower = AnaCVA - mat2(:,2); 
CVAUpper = AnaCVA + mat2(:,2); 


figure()
hold on 
plot( [0 50], [PD PD] , 'k-') 
plot(1:50, mat1(:,1), 'b.') 
plot(1:50, tauLower, 'b--') 
plot(1:50, tauUpper, 'b--') 
xlabel('Sample Size') 
ylabel('Probabilities') 
title('Estimates for the probability of default')
hold off 

figure() 
hold on 
plot( [0 50], [AnaCVA AnaCVA] , 'k-') 
plot(1:50, mat2(:,1), 'b.') 
plot(1:50, CVALower, 'b--') 
plot(1:50, CVAUpper, 'b--') 
xlabel('Sample Size') 
ylabel('Credit Value Adjustment') 
title('Estimates for the Credit Value Adjustment')
hold off 
