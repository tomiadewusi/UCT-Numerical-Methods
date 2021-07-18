clear 
clc
rng(0)

% Need to check this, I may be generating one extra random variate as 
% as the start and end points are fixed 
Z = randn(5,100);
% Z = randn(5,99); 
T = 5; 
t = 0:0.05:T;
S0 = 100; 
S100 = 250; 
mu = 0.15; 
sigma = 0.2; 
W0 = 0; 
%W100 = 1.3315;
W100 = (log(S100/S0) - (mu - 0.5*sigma^2)*T)/sigma; 

W = zeros(5,101);
S = zeros(5,101);  

for idx = 1:5
    %This generates the correct graph as shown in the textbook 
	W(idx,:) = [BrownianBridge(t,W0,W100,Z(idx,:)) W100];
    
    % W(idx,:) = [W0 BrownianBridge(t,W0,W100,Z(idx,:)) W100];
    
    % However if I do it this second way, I end up getting a 
    % a repeated value of zero in the beginning
	S(idx,:) = S0*exp((mu - 0.5*sigma^2)*t + sigma.*W(idx,:));
end 

figure()
subplot(1,2,1)
plot(t,W.')
subplot(1,2,2)
plot(t,S.')
