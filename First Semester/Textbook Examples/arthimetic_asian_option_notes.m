clear 
clc
rng(0)

% generating one asian option path 

T = 2; 
r = 0.085; 
sigma = 0.3;
K = 100; 
S0 = 100; 
t = linspace(0,2,1000);
deltat = t(2);
mat = zeros(50,1); 
for idx = 1:50
    N = 1000*idx; 
    Z = randn(length(t)-1,N);
    Si = [S0.*ones(1,N); S0*exp(cumsum((r-0.5*sigma.^2)*deltat...
        + sigma*sqrt(deltat).*Z))];

    option_value = mean(exp(-r*T)*max(  mean(Si)- K, 0)); 
    mat(idx,1) = option_value; 
end 
ave_option_value = mean(mat(11:end,1)); 
figure()
hold on
plot(1:50,mat)
plot([0 50],[ave_option_value ave_option_value])