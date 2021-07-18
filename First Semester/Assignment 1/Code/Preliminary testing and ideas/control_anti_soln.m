function value = control_anti_soln(sigma, K, eta, n)
% Fixed Input values
S0 = 100; L = 70; U = 130; r = 0.1; T = 1; N = 25; deltat = T/(N+1);

% n0  should be 10% of the total number of samples 

% Calculation of alpha
Z = randn(1000,N+1);
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(:,end );
Si = Si(:,1:end-1);
vanilla = exp(-r*T).*(max(eta.*(ST - K),0));
fade = vanilla.*mean(Si > L & Si < U,2);

Z = -Z; 
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(:,end );
Si = Si(:,1:end-1);
vanilla_m = exp(-r*T).*(max(eta.*(ST - K),0));
fade_m = vanilla_m.*mean(Si > L & Si < U,2); 

fade_ave = 0.5*(fade_m + fade); 
vanilla_ave = 0.5*(vanilla + vanilla_m); 

cov_fg = cov(fade_ave,vanilla_ave);
alpha = -cov_fg(1,2)/cov_fg(2,2);

% Generate standard normal variates and compute Antithetic estimates
Z = randn(n,N+1);
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(:,end );
Si = Si(:,1:end-1);
vanilla = exp(-r*T).*(max(eta.*(ST - K),0));
fade = vanilla.*mean(Si > L & Si < U,2);

Z = -Z; 
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(:,end );
Si = Si(:,1:end-1);
vanilla_m = exp(-r*T).*(max(eta.*(ST - K),0));
fade_m = vanilla_m.*mean(Si > L & Si < U,2); 

fade_ave = 0.5*(fade_m + fade); 
vanilla_ave = 0.5*(vanilla + vanilla_m); 

%Closed form solution
d1 = (log(S0/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
BSoption = eta*(S0*normcdf(eta*d1) - exp(-r*T)*K*normcdf(eta*d2));

% Final Value 
value = mean(fade_ave + alpha*(vanilla_ave - BSoption)); 
end 