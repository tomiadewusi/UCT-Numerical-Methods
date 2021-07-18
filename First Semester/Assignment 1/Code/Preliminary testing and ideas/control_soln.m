function value = control_soln(sigma, K, eta, n)
% Fixed Input values
S0 = 100; L = 70; U = 130; r = 0.1; T = 1; N = 25; deltat = T/(N+1);

% n0  should be 10% of the total number of samples 

%Closed form solution
d1 = (log(S0/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
BSoption = eta*(S0*normcdf(eta*d1) - exp(-r*T)*K*normcdf(eta*d2));

% Calculation of alpha
Z = randn(1000,N+1);
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(:,end );
Si = Si(:,1:end-1);
vanilla = exp(-r*T).*(max(eta.*(ST - K),0));
fade = vanilla.*mean(Si > L & Si < U,2);
cov_fg = cov(fade,vanilla);
alpha = -cov_fg(1,2)/cov_fg(2,2);

% Regenerate standard normal variates 
Z = randn(n,N+1);
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(:,end );
Si = Si(:,1:end-1);

% Crude Solution - vanilla
payoff_vanilla = exp(-r*T).*(max(eta.*(ST - K),0));

% Control Variates Solution
payoff_fade = payoff_vanilla.*mean(Si > L & Si < U,2);

vanilla_option = mean(payoff_vanilla);
fade_option = mean(payoff_fade);
value = fade_option + alpha*(vanilla_option - BSoption); 
end 