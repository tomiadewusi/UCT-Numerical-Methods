clear 
clc 

K = 100; 
A = 10; 
T = 1; 
r = 0.05; 
sigma = 0.2; 
L = 51; 

Sdom = linspace(50,150,L); 
tdom = linspace(0,T,L); 

[S,t] = meshgrid(Sdom,tdom); 

C = zeros(size(t));

d = (log(S(t<T)./K)+(r-0.5.*sigma.^2).*(T-t(t<T)))./(sigma.*sqrt(T-t(t<T))); 
N = 0.5.*(1+erf(d./sqrt(2))); 


C(t<T) = A.*exp(-r.*(T-t(t<T))).*N; 
C(t==T) = 0.5.*A.*(1 + sign(S(t==T)-K)); 

figure() 
surf(S,t,C)
title('Surface for a cash-or-nothing Call option'); 
xlabel('Strike Price')
ylabel('Time')
zlabel('Option value')


