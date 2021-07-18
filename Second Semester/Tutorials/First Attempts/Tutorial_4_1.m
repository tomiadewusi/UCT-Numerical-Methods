clear 
clc  
rng(0) 

S0 = 40; 
sigma = 0.6; 
alpha = 0.85; 
r = 0.06; 
deltaT = 0.001; 
deltaK = 0.001; 

kappa = @(T) (2*r)./( (sigma.^2).*(1-alpha).*(exp(2*r*(1-alpha).*T)-1) );
x = @(T) kappa(T).*S0.^(2*(1-alpha)).*(exp(2*r*(1-alpha).*T)); 
y = @(T,K) kappa(T)*K.^(2.*(1-alpha));
z = 2 + 1/(1-alpha);

CallExact = @(T,K) S0.*(1 - ncx2cdf(y(T,K),z,x(T)))...
    - K.*exp(-r*T).*ncx2cdf(x(T),z-2,y(T,K)); 

PutExact = @(T,K) -S0.*ncx2cdf(y(T,K),z,x(T))...
    + K.*exp(-r*T).*(1-ncx2cdf(x(T),z-2,y(T,K))); 

dPdT = @(T,K) (PutExact(T+deltaT,K)-PutExact(T-deltaT,K))./...
    (2*deltaT) ; 

dPdK = @(T,K) (PutExact(T,K+deltaK)-PutExact(T,K-deltaK))./...
    (2*deltaK) ; 

d2PdK2 = @(T,K) (   PutExact(T,K+2*deltaK)...
    -2.*PutExact(T,K) + PutExact(T,K-2*deltaK)  )./...
    (4*deltaK.^2) ; 

dCdT = @(T,K) (CallExact(T+deltaT,K)-CallExact(T-deltaT,K))./...
    (2*deltaT) ; 

dCdK = @(T,K) (CallExact(T,K+deltaK)-CallExact(T,K-deltaK))./...
    (2*deltaK) ; 

d2CdK2 = @(T,K) (   CallExact(T,K+2*deltaK)...
    -2.*CallExact(T,K) + CallExact(T,K-2*deltaK)  )./...
    (4*deltaK.^2) ; 

DupirePut = @(T,K) (sqrt(2)./K).*(sqrt((dPdT(T,K)+r.*K.*dPdK(T,K))...
    ./d2PdK2(T,K))); 
DupireCall = @(T,K) (sqrt(2)./K).*(sqrt((dCdT(T,K)+r.*K.*dCdK(T,K))...
    ./d2CdK2(T,K))); 
CEVSigma = @(S) sigma.*S.^(alpha-1); 

T = 2; 
figure()
hold on 
DupireSoln = [DupireCall(T,1:39),DupirePut(T,40:80)]; 
plot(1:80,DupireSoln,'r.')
plot(1:80,CEVSigma(1:80)) 
hold off 

figure()
plot(1:80,(DupireSoln - CEVSigma(1:80)))
ylim([-5e-5 5e-5])
