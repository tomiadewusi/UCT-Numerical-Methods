clear 
clc

%Step 1 - assign variables 

K = 100; 
A = 10; 
T = 1; 
r = 0.05; 
sigma = 0.2; 
L = 51; 

%Step 2 -  Sdom

sdom = linspace(50, 150, L); 
tdom = linspace(0, T, L); 

%Step 3 - meshgrid 

[S, t] = meshgrid(sdom, tdom); 

[number] = size(S); 

%Step 4 
%Even though when you index d and n using the logical, it comes out as a
%vector, when you go into into C it sorts it self out 
C = zeros(size(t)); 

d = (  log(S(t<T)/K) + (r - 0.5.*sigma.^2 ).*(T-t(t<T))  )./(  sigma.*sqrt(T-t(t<T)) ); 
N = 0.5.*(1 + erf(d./sqrt(2))); 
C(t==T) = 0.5.*A.*(1 + sign(S(t==T)-K)); 
C(t<T) = A.*exp(-r.*(T-t(t<T))).*N;



figure()

surf(S,t,C)
title('Cash or nothing Call Option')
xlabel('Strike Price')
ylabel('Time')
zlabel('Call option value')






