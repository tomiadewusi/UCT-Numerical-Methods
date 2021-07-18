function value = anti_soln(sigma, K, eta, n)
% Fixed input values 
S0 = 100; L = 70; U = 130; r = 0.1; T = 1; N = 25; deltat = T/(N+1);

Z = randn(n,N+1);
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(:,end );
Si = Si(:,1:end-1);
payoff_fade = exp(-r*T).*(max(eta.*(ST - K),0)).*mean(Si > L & Si < U,2);

Z = -Z; 
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(:,end );
Si = Si(:,1:end-1);
payoff_fade_m = exp(-r*T).*(max(eta.*(ST - K),0)).*mean(Si > L & Si < U,2);

value = 0.5*(mean(payoff_fade + payoff_fade_m)); 
end 

