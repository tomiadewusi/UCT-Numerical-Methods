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

% Generating the closed sequence 
Z2 = vanderCorput(50001,2);
Z3 = vanderCorput(50001,3); 
Z5 = vanderCorput(50001,5); 


mat1 = zeros(50,1);
mat2 = zeros(50,1);
for idx = 1:50
    % Crude Monte-Carlo
    n = 1000*idx;
    Z0 = randn(n,N);
    index = ((1:n)./(n+1)).'; 
    ham = [index Z2(2:n+1) Z3(2:n+1) Z5(2:n+1)];
    Z = norminv(ham); 
    
    inside = (r - 0.5*sigma^2)*h + sigma*sqrt(h).*Z;
    Si = [S0.*ones(n,1) S0*exp(cumsum(inside,2))];
    sample = max(diff(Si,1,2),0)*exp(-r*h.*(1:N).');
    option = mean(sample);
    mat1(idx,:) = option;
    
    % Monte-Carlo with end point stratification  
    % The Brownian Bridge
    Z = norminv(index); 
    ham1 = [Z3(2:n+1) Z2(2:n+1) Z5(2:n+1)]; 
    Z1 = norminv(ham1); % generating all at the same time 
    ST = S0*exp((r - 0.5*sigma^2)*T + sigma*sqrt(T).*Z); 
    W0 = zeros(n,1); 
    W100 = ((log(ST./S0) - (r - 0.5*sigma^2)*T)/sigma); 
    W = [W0 BrownianBridgeV(t,W0,W100,Z1 ) W100];  
    
    Si = S0*exp((r - 0.5*sigma^2)*t + sigma.*W);
    fx = max(diff(Si,1,2),0)*exp(-r*h.*(1:N).');
    option = mean(fx); 
    mat2(idx,:) = option;
end

figure()
hold on 
plot([0 50], [c c] , 'k-')
plot(mat1(:,1), 'b.')
plot(mat2(:,1), 'r.')
hold off 
