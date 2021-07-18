clear 
clc 

f = @(x) exp(-((x.^2)./2))./(sqrt(2.*pi)); 
a = -3; 
b = 3; 
h = 0.1; 

format long 
%you need to very careful with variable names 
c = composite(f,a,b,h); 
m = midpoint(f,a,b,h); 
inte = integral(f,a,b); 
fo = forward(f,a,b,h); 
cen = central(f,a,b,h); 
points = a:h:b; 
fun = f(points); 
%you need to read the documentation carefully to avoid mistakes 
grad = gradient(f(points),h); 

figure() 
Y = [fun' grad' cen' fo']; 
X = repmat(points', 1, 4); 
plot(X,Y)
legend('Function', 'Gradient', 'Central', 'Forward') 


