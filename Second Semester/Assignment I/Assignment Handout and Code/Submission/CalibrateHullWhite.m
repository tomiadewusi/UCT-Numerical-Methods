function x = CalibrateHullWhite
% ADWOLA002 - Tomi Adewusi 
%--------------------------------------------------------------------------
%% Data
%--------------------------------------------------------------------------
N=20;
Dates=0:N+1; % Times T_0 to T_{N+1}
% Bond prices for T_0 to T_{N+1}
Bstar=[1.000000000000000 0.931733068514009 0.867336864194859 0.807180507606257 0.751390186399840 ...
    0.699924481765438 0.652631479135652 0.609290924427683 0.569644566898682 0.533417438303370 ...
    0.500332340416184 0.470119346466237 0.442521711052162 0.417299242183633 0.394229917259718 ...
    0.373110314043405 0.353755267713142 0.335997045644106 0.319684243642614 0.304680543400063 ...
    0.290863424934260 0.278122895071373];
T_i = Dates; 
r0 = 0.07; 
% Caplet prices with strike K=0.05 and maturities T_1 to T_{N-1}
Caplets=[0.022922737518613 0.023535825778635 0.023194289857391 0.022294316394003 0.021086407641792 ...
    0.019728461392470 0.018320943715610 0.016927232296721 0.015586172322853 0.014320233052482 ...
    0.013140971196705 0.012052749925964 0.011055287463489 0.010145403965974 0.009318213062292 ...
    0.008567926766803 0.007888391002214 0.007273433858145 0.006717084341419];
F = griddedInterpolant(T_i,-log(Bstar)); 
% Ti_trunc = Dates(1:end-1);
Ti_trunc = 1:20; 
delta = 0.00000000001; 
fwd_rates = [r0, (F(Ti_trunc + delta) - F(Ti_trunc - delta))/(2*delta)];
fun =@(x) CaliObjective(Bstar,fwd_rates,Caplets,x); 
options = optimset('TolX',1e-50,...
                   'TolFun',1e-50,...
                   'Display','off');     
x  = fminsearch(fun,[0.1 0.1],options);
end

function error = CaliObjective(Bstar,fwd_rates,Caplets,x)
alpha = x(1); 
sigma = x(2); 
K = 0.05; 
r0 = 0.07; 
Ti = 1:19; 
% These are row vectors 
output = CapletAnalytical(Bstar,fwd_rates,Ti,K,r0,alpha,sigma)-Caplets; 
error = output*output.';  
end

function Price = CapletAnalytical(Bstar,fstar,Ti,K,r0,alpha,sigma)
Price = (1+K).*PutPrice(Bstar,fstar,Ti, Ti+1, 1/(1+K),r0,alpha,sigma); 
end

function Price = PutPrice(Bstar,fstar, t1,t2,kappa,r0,alpha,sigma)

A = @(t1,t2)  (1./alpha).*(1-exp(-alpha.*(t2-t1))); 

C = @(t1,t2) log(Bstar(t2+1)./Bstar(t1+1)) ...
       + fstar(t1+1).*A(t1,t2) ...
       - ((sigma.^2)/(4*alpha)).*(1-exp(-2*alpha.*t1)).*A(t1,t2).^2 ; 
   
B = @(t1,t2,r1) exp( - A(t1,t2).*r1 + C(t1,t2) ); 

sigma_ave = sigma.*sqrt((1-exp(-2*alpha.*t1))./(2*alpha) ).*A(t1,t2); 

d1 = (log(B(0,t2,r0)./(B(0,t1,r0).*kappa)) + 0.5.*sigma_ave.^2)./sigma_ave; 
d2 = d1 - sigma_ave; 
Price = kappa.*B(0,t1,r0).*normcdf(-d2) - B(0,t2,r0).*normcdf(-d1);
end 

