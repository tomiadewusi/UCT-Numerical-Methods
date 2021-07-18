clear 
clc 
rng(0)
T = 2; 
r = 0.085; 
sigma = 0.3;
K = 100; 
S0 = 100; 

% closed form solution 
sigmaG = sigma/sqrt(3); 
b = 0.5*(r - 0.5*sigmaG.^2); 
d1 = ( log(S0/K) + (b + 0.5*sigmaG.^2)*T )/( sigmaG*sqrt(T)); 
d2 = d1 - sigmaG*sqrt(T); 

callG = S0*exp((b-r)*T)*normcdf(d1) - K*exp(-r*T)*normcdf(d2); 

t = linspace(0,2,1000);
deltat = t(2);
mat = zeros(50,1); 

for idx = 1:50
    N = 1000*idx; 
    Z = randn(length(t)-1,N);
    Si = [S0.*ones(1,N); S0*exp(cumsum((r-0.5*sigma.^2)*deltat...
        + sigma*sqrt(deltat).*Z))];

    option_value = mean(exp(-r*T)*max( geomean(Si)- K, 0)); 
    mat(idx,1) = option_value; 
end 

figure()
hold on
plot(1:50,mat,'k.')
plot([0 50],[callG callG])