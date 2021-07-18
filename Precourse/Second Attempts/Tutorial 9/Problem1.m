clear 
clc 

h = 0.025; 
t = 0:h:10;
len = length(t); 
W = zeros(size(t)); 
w_0 = 1; 
W(1) = w_0; 
w_old = w_0; 
f = @(t,y) 15.*t.^2 - 15.*y + 2.*t; 


for indexi = 2:len
    w_new = w_old + h.*f(t(indexi - 1) + (h./2),...
        w_old + (h./2).*f(t(indexi - 1), w_old)); 
    w_old = w_new; 
    W(indexi) = w_new; 
end 

[t45,W_ode45] = ode45(f, t, w_0); 
[t15s,W_ode15s] = ode15s(f, t, w_0); 

y = exp(-15.*t) + t.^2; 

figure()
subplot(2,3,1)
plot(t, W)
title('Midpoint')

subplot(2,3,4)
plot(t, W - y) 
title('Midpoint Error')

subplot(2,3,2)
plot(t45, W_ode45)
title('Ode45') 

subplot(2,3,5)
plot(t45, W_ode45 - y')
title('Ode45 Error') 

subplot(2,3,3)
plot(t15s, W_ode15s)
title('Ode15s')

subplot(2,3,6)
plot(t15s, W_ode15s - y')
title('Ode15s Error')

