clear
clc
rng(0)

S0 = 50;
K = 55;
sigma = 0.25;
r = 0.1;
T = 2;
t = 0;
n = 100000; 
delta = 0.2; % Recovery Rate

BHat = 75;
Vt = 100;
spread = 0.025;
Dt = (BHat*exp(-r*(T-t)))/(Vt);

fun = @(sigmaHat) spread - (-1/(T-t))*...
    (log(    (1/Dt)...
    *normcdf(-((log(Vt/BHat)+(r+0.5*sigmaHat.^2)*(T-t))...
    /(sigmaHat*sqrt(T-t))))...
    +normcdf((log(Vt/BHat)+(r-0.5*sigmaHat.^2)*(T-t))...
    /(sigmaHat*sqrt(T-t)))    )) ;

sigmaHat = fzero(fun,sigma);
mat1 = zeros(21,2);
mat2 = zeros(21,2);
for ix = 1:21
    rho = -1 + (ix-1)*0.1;
    z = randn(2,n);
    z1 = z(1,:);
    z2 = z(2,:);
    zv = rho.*z1 + sqrt(1 - rho.^2).*z2;
    ST = S0*exp((r - 0.5*sigma.^2)*(T-t) + sigma*sqrt(T-t).*z1);
    VT = Vt*exp((r - 0.5*sigmaHat.^2)*(T-t) + sigmaHat*sqrt(T-t).*zv);
    Indicator = VT < BHat;  % Default or nah
    Vmean = mean(Indicator); % Probability of Default
    Vconf = 3*std(Indicator)/sqrt(n);
    mat1(ix,:) = [Vmean, Vconf]; 
    CVA = exp(-r*(T-t))*(1-delta).*max(ST - K,0).*Indicator; 
    mat2(ix,:) = [mean(CVA), 3*std(CVA)/sqrt(n)] ; 
end
%Analytical Solution 
PD = normcdf(-( (log(Vt/BHat)+(r-0.5*sigmaHat.^2)*(T-t))...
    /(sigmaHat*sqrt(T-t)) ));
d1 = ((log(S0/K)+(r+0.5*sigma.^2)*(T-t))...
    /(sigma*sqrt(T-t)));
d2 = ((log(S0/K)+(r-0.5*sigma.^2)*(T-t))...
    /(sigma*sqrt(T-t)));

X0 = S0*normcdf(d1) - K*exp(-r*(T-t))*normcdf(d2); % Assuming that default
% can only occur at a the terminal time
AnaCVA = (1-delta)*X0*PD; 

Vlower = PD - mat1(:,2) ; 
Vupper = PD + mat1(:,2) ; 

CVAlower = mat2(:,1) - mat2(:,2) ; 
CVAupper = mat2(:,1) + mat2(:,2) ; 


figure()
hold on 
plot( (-1:0.1:1), mat1(:,1), 'b.') 
plot( (-1:0.1:1)', Vupper, 'b--') 
plot( (-1:0.1:1)', Vlower, 'b--') 
plot( [-1 1], [PD PD], 'k-')
title('Estimates for the Probablity of default')
xlabel('Correlation')
ylabel('Probabilites')
hold off 

figure() 
hold on 
plot( (-1:0.1:1), mat2(:,1), 'b.') 
plot( (-1:0.1:1)', CVAupper, 'b--') 
plot( (-1:0.1:1)', CVAlower, 'b--') 
plot( [-1 1], [AnaCVA AnaCVA], 'k-')
xlabel('Correlation')
ylabel('Credit Value Adjustment')
hold off 