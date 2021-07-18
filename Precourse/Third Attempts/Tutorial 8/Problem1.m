clear 
clc 

f = @(x) exp(-(x.^2)./2)./(sqrt(2.*pi)); 
a = -3; 
b = 3; 
h = 0.1; 
t = a:h:b; 

v1 = midpoint(f,a,b,h);
v2 = composite(f,a,b,h);
v3 = integral(f,a,b); 

v4 = forw(f,a,b,h); 
v5 = central(f,a,b,h); 
v6 = gradient(f(t),h); 

figure()
plot([t.' t.' t.' t.'], [f(t).' v6.' v4.' v5.'])
