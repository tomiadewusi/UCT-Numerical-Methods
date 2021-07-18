clear 
clc
% rng(1, 'v4') 

%This parameter changes the number of brownian motion paths 
k = 100;

mu = 0.3; 
sigma = 0.1; 
X_0 = 10; 
delta_t = 0.01; 
t = 0:delta_t:5;
[~,t_length] = size(t); 
numbers = randn(t_length-1, k);

tic 
%by default, if provided with a matrix then cumprod works along the columns
temp = X_0.*cumprod( exp( (mu - 0.5.*sigma.^2).*delta_t +...
    sigma.*sqrt(delta_t).*numbers )); 
X = [repmat(X_0,1, k) ; temp ]; 
T = repmat(t', 1, k); 

%When plotting two matrices together, matlab takes in pairs of columns 
figure()
hold on 
plot(T,X)
plot(t, X_0.*exp(mu.*t),'k-')
%axis tight
hold off 
temp = X_0.*exp(mu.*t); 
count = sum(X(end,:) > temp(end)); 
%I am checking the number of paths that are above the expected value 
count/k 
toc 