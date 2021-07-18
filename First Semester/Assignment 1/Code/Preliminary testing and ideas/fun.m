
function value = fun(Z)

K = 80; sigma = 0.2677; eta = 1;

S0 = 100; L = 70; U = 130; r = 0.1; T = 1; N = 25; deltat = T/(N+1);
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(end);
Si = Si(1:end-1);
value = -(log(mean(exp(-r*T).*(max(eta.*(ST - K),0)).*mean(Si > L & Si < U,2)))...
    -0.5.*Z*Z.');
end