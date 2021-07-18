clear 
clc
rng(0) 

mu = 0.5; 
sigma = 0.4; 
S0 = 1; 
alpha = 1.2; 
t = linspace(0,5,1000);  
deltat = t(2) - t(1); 
Z = randn(length(t)-1, 1); 

S_euler = zeros(length(t),1); 
S_euler(1) = S0;
S_mils = zeros(length(t),1); 
S_mils(1) = S0; 

for idx = 2:length(t)
	S_euler(idx) = S_euler(idx-1)+ mu*S_euler(idx-1).*deltat...
		+ (sigma*S_euler(idx-1)^(alpha*0.5))*Z(idx-1)*sqrt(deltat); 
end 

for idx = 2:length(t)
	S_mils(idx) = S_mils(idx-1)+ mu*S_mils(idx-1).*deltat...
		+ (sigma*S_mils(idx-1)^(alpha*0.5))*Z(idx-1)*sqrt(deltat)...
        + (alpha*0.250*sigma^2.*(S_mils(idx-1)^(alpha-1)).*(Z(idx-1)^2 - 1)...
        .*deltat); 
end 

figure()
subplot(1,3,1)
plot(t,S_euler)
title('The Euler-Maruyama Approximation')

subplot(1,3,2)
plot(t,S_mils)
title('The Milstein Approximation')

subplot(1,3,3)
plot(t,S_euler - S_mils) 
title('The Difference: Euler - Mils ')



