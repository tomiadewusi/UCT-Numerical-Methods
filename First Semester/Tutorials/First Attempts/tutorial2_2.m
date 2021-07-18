clear 
clc 
rng(0) 


%this is used to generate from the importance function 
%Using the probability integral transform, we can get
%realisations from an exponential with rate 1. This  
Y = @(u) -1*log(u); 
iterations = 10000; 
X = zeros(iterations ,1); 
counter = 0; 

while  counter ~= iterations 
	u = rand(2,1);
	if u(2) > exp( -0.5*( (Y(u(1)) -1)^2  ) )
		continue
	else 
		counter = counter + 1;
		%this should only be simulated if you enter here 
		v = rand(1,1); 
		if v < 0.5
			X(counter) = -Y(u(1));
		else 
			X(counter) = Y(u(1));
		end 
	end 
end 
 


figure(1)
subplot(1,2,1)
hist(X) 
subplot(1,2,2)
hist(randn(iterations,1))

#Extra Plots from the Notes
figure(2) 
subplot(1,2,2)
t = linspace(0,6,1000);
plot(t,exp( -0.5.*( (t -1).^2  ) ) ,"k-")
ylabel("Conditional probability of accepting an exponential generated sample")
title("Conditional probability of acceptance")

subplot(1,2,1)
hold on 
plot(t,2./sqrt(2.*pi).*exp(-0.5.*t.^2),'k-')
plot(t, exp(-t).*sqrt(2*exp(1)/pi),'k--')