clear 
clc 
% rng(0) 

% Deriving and comparing the Pathwise estimates for all the Greeks 
% https://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model
% Jump to the section on the greeks for all the formulas

%--------------------------------------------------------------------------

% Fixed Inputs 
S0 = 100; 
K = 80; 
r = 0.1; 
sigma = 0.1; 
T = 2; % This is the maturity time 
t = 0.;  % This is the current time 
tau = T-t;  % This is the time to maturity 
eta = 1; 

%--------------------------------------------------------------------------
% This is the analytical solution 
d1 = (log(S0/K) + (r + 0.5*sigma^2).*tau)/(sigma*sqrt(tau)); 
d2 = d1 - sigma*sqrt(tau); 
c = eta*(S0*normcdf(eta*d1) - normcdf(eta*d2)*K*exp(-r*tau)); 

% Derivative with respect to the underlying S0
Delta = eta*normcdf(eta*d1); 

% Second derivative w.r.t. the underlying or Derivative w.r.t. to Delta 
Gamma = normpdf(d1)/(S0*sigma*sqrt(tau)); 

% Derivative w.r.t. to the volatility - sigma 
Vega = S0*normpdf(d1)*sqrt(tau); 

% Gamma and Vega are the same for calls and puts -> Put-Call-Parity

% Derivative w.r.t. time to maturity - T-t
Theta = -((S0*normpdf(d1)*sigma)/(2*sqrt(tau)))...
    -eta*r*K*exp(-r*tau)*normcdf(eta*d2); 

% Derivative w.r.t. the risk free rate
Rho = eta*K*tau*exp(-r*tau)*normcdf(eta*d2); 

%--------------------------------------------------------------------------
% Monte Carlo Simulation 
% Forward Difference Methods 
n = 10; 
Z1 = randn(n,1); 
Z2 = randn(n,1); 
S0Bump = 0.1; 
f1 = @(Z) exp(-r*tau)*max(S0*exp((r-0.5*sigma.^2)*tau...
    + sigma*sqrt(tau).*Z) - K,0) ;  
f2 = @(Z) exp(-r*tau)*max((S0+S0Bump)*exp((r-0.5*sigma.^2)*tau...
    + sigma*sqrt(tau).*Z) - K,0);  

% DeltaFD = mean(f2(Z2) - f1(Z1))/S0Bump; 
% DeltaFD_CN = mean(f2(Z2) - f1(Z2))/S0Bump; 
% Delta; 

% Apply the same principle to other greeks or multiple derivatives

%--------------------------------------------------------------------------
% Pathwise Derivative Estimates
% Delta
ST = @(Z) S0*exp((r-0.5*sigma.^2)*tau + sigma*sqrt(tau).*Z)  ;  
DeltaP = mean(exp(-r*tau).*(ST(Z1) > K)...
    .*ST(Z1)./S0); 
Delta; 
GammaP = mean(exp(-r*tau).*(ST(Z1) > K)...
    .*ST(Z1)./(S0.^2))
Gamma
%--------------------------------------------------------------------------
% Vega 

% VegaP = mean(exp(-r*tau).*(ST(Z1) > K)...
%     .*ST(Z1).* ((log(ST(Z1)./S0) - (r+0.5*sigma.^2)*tau)/sigma))

% Bump = 0.01;
% VegaP1 = mean(exp(-r*tau).*(ST(Z1) > K)...
%     .*ST(Z1).*(-sigma*tau + sqrt(tau).*Z1)) 
% 
% g2 = @(Z) exp(-r*tau)*max(S0*exp((r-0.5*(sigma + Bump).^2)*tau...
%     + (sigma + Bump)*sqrt(tau).*Z) - K,0) ;  
% 
% g1 = @(Z) exp(-r*tau)*max(S0*exp((r-0.5*sigma.^2)*tau...
%     + sigma*sqrt(tau).*Z) - K,0) ; 
% 
% VegaFD_CN = mean(g2(Z2) - g1(Z2))/Bump
% Vega 


