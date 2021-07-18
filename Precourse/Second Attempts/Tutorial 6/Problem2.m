clear 
clc

K = 100; 
A = 10; 
T = 1; 
r = 0.05; 
sigma = 0.2;
L = 51; 

Sdom = linspace(50,150,L); 
tdom = linspace(0,T,51); 

[S,t] = meshgrid(Sdom,tdom); 
C = zeros(size(t)); 

d = ( log(S./K) + (r - 0.5.*sigma.^2).*(T-t) )./(sigma.*sqrt(T - t)); 

N = 0.5.*(1 + erf(d./sqrt(2))); 

C(t<T) = A.*exp(-r.*(T - t(t<T))).*N(t<T); 

C(t == T) = 0.5.*(1 + sign(S(t == T) - K)); 

figure()
surf(S,t,C)
ylabel("Time")
xlabel("Strike Price")
zlabel("Option Value")
title("Surface for a cash or nothing call option")


