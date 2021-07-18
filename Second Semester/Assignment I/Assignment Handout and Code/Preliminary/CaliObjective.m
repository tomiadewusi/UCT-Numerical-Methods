function error = CaliObjective(Bstar,fwd_rates,Caplets,x)
alpha = x(1); 
sigma = x(2); 
K = 0.05; 
r0 = 0.07; 
Ti = 1:19; 
% These are row vectors 
output = CapletAnalytical(Bstar,fwd_rates,Ti,K,r0,alpha,sigma); 
error = (output - Caplets)*(output - Caplets).'; 
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

