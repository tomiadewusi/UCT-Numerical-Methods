function value = lhs_soln(sigma, K, eta, n)
% Fixed Input values
S0 = 100; L = 70; U = 130; r = 0.1; T = 1; N = 25; deltat = T/(N+1);
Z = lhsnorm(zeros(1,N+1),eye(N+1),n); 
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(:,end );
Si = Si(:,1:end-1);
value = mean(exp(-r*T).*(max(eta.*(ST - K),0)).*mean(Si > L & Si < U,2));
end 