clear
clc

N=20;
Dates=0:N+1; % Times T_0 to T_{N+1}
% Bond prices for T_0 to T_{N+1}
Bstar=[1.000000000000000 0.931733068514009 0.867336864194859 0.807180507606257 0.751390186399840 ...
0.699924481765438 0.652631479135652 0.609290924427683 0.569644566898682 0.533417438303370 ...
0.500332340416184 0.470119346466237 0.442521711052162 0.417299242183633 0.394229917259718 ...
0.373110314043405 0.353755267713142 0.335997045644106 0.319684243642614 0.304680543400063 ...
0.290863424934260 0.278122895071373];

% Caplet prices with strike K=0.05 and maturities T_1 to T_{N-1}
Caplets=[0.022922737518613 0.023535825778635 0.023194289857391 0.022294316394003 0.021086407641792 ...
0.019728461392470 0.018320943715610 0.016927232296721 0.015586172322853 0.014320233052482 ...
0.013140971196705 0.012052749925964 0.011055287463489 0.010145403965974 0.009318213062292 ...
0.008567926766803 0.007888391002214 0.007273433858145 0.006717084341419];

[x,fval,exitflag,output] = CalibrateHullWhite(Bstar,Caplets,Dates); 
exitflag
r0 = 0.07; 
t1 = zeros(1,22); 
t2 = Dates; 
F = griddedInterpolant(t2,-log(Bstar)); 
delta = 0.00000000001; 
Ti_trunc = Dates(2:end-1); 
fstar = [r0, (F(Ti_trunc + delta) - F(Ti_trunc - delta))/(2*delta)]; 

CaliObjective(Bstar,fstar,Caplets,x)
alpha = x(1)
sigma = x(2)

A = @(t1,t2)  (1./alpha).*(1-exp(-alpha.*(t2-t1))); 

C = @(t1,t2) log(Bstar(t2+1)./Bstar(t1+1)) ...
       + fstar(t1+1).*A(t1,t2) ...
       - ((sigma.^2)/(4*alpha)).*(1-exp(-2*alpha.*t1)).*A(t1,t2).^2 ; 
   
B = @(t1,t2,r1) exp( - A(t1,t2).*r1 + C(t1,t2) );
% This shows that under those values of alpha and sigma, we are able to
% reproduce the bond curve to a high degree of accuracy. 
% [Bstar.' B(t1,t2,r0).']
figure()
subplot(1,2,1)
plot(t2,(Bstar-B(t1,t2,r0)) )
subplot(1,2,2)
plot(t2,log10(abs(Bstar-B(t1,t2,r0))))
