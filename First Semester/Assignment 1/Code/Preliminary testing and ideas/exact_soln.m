function value = exact_soln(sigma, K, eta)
% Fixed Input values
S0 = 100; L = 70; U = 130; r = 0.1; T = 1; N = 25; deltat = T/(N+1);
t = deltat.*(1:25);

d1 = (log(S0/K) + (r+0.5*sigma^2)*T )/(sigma*sqrt(T)); 
d2 = d1 - sigma*sqrt(T); 

d3 = (log(S0/L) + (r+0.5*sigma^2).*t )./(sigma.*sqrt(t)); 
d4 = d3 - sigma.*sqrt(t); 

d5 = (log(S0/U) + (r+0.5*sigma^2).*t )./(sigma.*sqrt(t)); 
d6 = d5 - sigma.*sqrt(t); 

M = @(x,y,rho) mvncdf([x y],0,[1 rho; rho 1]);
p = -eta.*sqrt(t)./sqrt(T); 

tot = 0; 
for ix = 1:N 
    tot = tot + eta*(S0*(M(-d5(ix),eta*d1,p(ix))-M(-d3(ix), eta*d1,p(ix)))...
    -K*exp(-r*T)*(M(-d6(ix),eta*d2,p(ix))-M(-d4(ix), eta*d2,p(ix)))); 
end 
value = tot/N; 
end