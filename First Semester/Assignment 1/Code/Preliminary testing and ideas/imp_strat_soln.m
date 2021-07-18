function value = imp_strat_soln(K,sigma,eta,n)
S0 = 100; L = 70; U = 130; r = 0.1; T = 1; N = 25; deltat = T/(N+1);
f = @(Z) fun(Z,K,sigma,eta,S0,L,U,r,T,deltat);
Z0 = ones(1,N+1);
v = fminsearch(f,Z0,options);
v = v/norm(v);
v = v.'
K1 = n;
final = zeros(K1,N+1); 
for idx = 1:K1
    U1 = rand(1);
    V = (idx - 1 + U1)/K1;
    X = norminv(V);
    Z = randn(N+1,1);
    I = eye(N+1,N+1);
    final(idx,:) = I*v.*X + (I - I*(v*v.')*I)*Z; 
end
Z = final;
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(:,end ); 
Si = Si(:,1:end-1);
value = mean(exp(-r*T).*(max(eta.*(ST - K),0)).*mean(Si > L & Si < U,2));
end 

function value = fun(Z,K,sigma,eta,S0,L,U,r,T,deltat)
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(end);
Si = Si(1:end-1);
value = -(log(mean(exp(-r*T).*(max(eta.*(ST - K),0)).*mean(Si > L & Si < U,2)))...
    -0.5.*Z*Z.');
end