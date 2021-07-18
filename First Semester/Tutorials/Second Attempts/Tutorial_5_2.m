clear 
clc 
rng(0) 
% Need to double check this with the boys
N = 1000; 
T = 5; 
t = linspace(0,T,N+1); 
S0 = 1; 
mu = 0.5; 
sigma = 0.4; 
alpha = 1.2; 

Sem = zeros(N+1,1); 
Sem(1) = S0; 
Smil = zeros(N+1,1);
Smil(1) = S0; 
Z = randn(N,1); 

for ix = 2:N+1
    Sem(ix) = Sem(ix-1) + mu*Sem(ix-1)*t(2) ...
        + sigma*Z(ix-1)*sqrt(t(2))*Sem(ix-1)^(alpha/2);
    
    Smil(ix) = Smil(ix-1) + mu*Smil(ix-1)*t(2) ...
        + sigma*Z(ix-1)*sqrt(t(2))*Smil(ix-1)^(alpha/2)...
        + ((alpha/4)*(sigma.^2)*(Smil(ix-1)^(alpha-1)))...
        *(Z(ix-1).^2*t(2) - t(2)) ; 
end 

figure() 
subplot(1,3,1) 
plot(t,Sem) 
subplot(1,3,2) 
plot(t,Smil)
subplot(1,3,3) 
plot(t, Sem - Smil) 