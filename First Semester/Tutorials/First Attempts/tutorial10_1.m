clear
clc
rng(0)

% Inputs
S0 = 100;
sigma = 0.3;
r = 0.1;
T = 1;

N = 4;
h = T/N;
d = 50;
t = 0:h:T; 
% Analytical Solution
d1 = ((r + 0.5*sigma^2)*h)/(sigma*sqrt(h));
d2 = d1 - sigma*sqrt(h);
c = N*S0*( normcdf(d1) - exp(-r*h)*normcdf(d2) );

mat1 = zeros(50,3);
mat2 = zeros(50,3);
for idx = 1:50
    % Crude Monte-Carlo
    n = 1000*idx;
    Z = randn(n,N);
    inside = (r - 0.5*sigma^2)*h + sigma*sqrt(h).*Z;
    Si = [S0.*ones(n,1) S0*exp(cumsum(inside,2))];
    sample = max(diff(Si,1,2),0)*exp(-r*h.*(1:N).'); 
    % Matrix multiplication gives you the sum  
    option = mean(sample);
    confidence = 3*std(sample)/sqrt(n);
    lower = c - confidence;
    upper = c + confidence;
    mat1(idx,:) = [option lower upper];
    % Monte-Carlo with end point stratification 
    % NB! The means and the variance are calculated per strata
    nums = (((1:d) - 1)/d) + (1/d).*rand(n/d,d);
    Z = norminv(nums); % Used for the Wiener process terminal values 
    % The Brownian Bridge
    Z = reshape(Z,[],1); 
    Z1 = randn(n,N-1); % generating all at the same time 
    ST = S0*exp((r - 0.5*sigma^2)*T + sigma*sqrt(T).*Z); 
    W0 = zeros(n,1); 
    W100 = ((log(ST./S0) - (r - 0.5*sigma^2)*T)/sigma); 
    W = [W0 BrownianBridgeV(t,W0,W100,Z1 ) W100];  
    
    Si = S0*exp((r - 0.5*sigma^2)*t + sigma.*W);
    fx = max(diff(Si,1,2),0)*exp(-r*h.*(1:N).');
    fx = reshape(fx, n/d,d); % so that suitable for Stratication 
    option = (1/d)*sum(mean(fx));
    option_var = (1/d)*sum(var(fx))/n; 
    confidence = 3*sqrt(option_var);
    lower = c - confidence;
    upper = c + confidence;
    mat2(idx,:) = [option lower upper];
end

figure()
hold on 
plot([0 50], [c c] , 'k-')
plot(mat1(:,1), 'b.')
plot(mat1(:,2:3), 'b--')
plot(mat2(:,1), 'r.')
plot(mat2(:,2:3), 'r--')
hold off 