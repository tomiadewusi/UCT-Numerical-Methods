function vector = HestCharST(u,S0,nu0,theta,kappa,sigma,rho,T,r)
%HestCharST Calculates the value of the Heston Characteristic Function 

a = 0.5*sigma.^2; % scalar
b = kappa - rho*sigma*1i.*u; % vector  
c = -0.5.*(u.^2 + 1i.*u); % vector
d = sqrt(b.^2 - 4*a.*c); % vector 
xmin = (b-d)./(2.*a); % vector 
xplu = (b+d)./(2.*a); % vector 
g = xmin./xplu; % vector 
D = xmin.*(1 - exp(-T.*d))./(1 - g.*exp(-T*d)); 
C = r*T*1i.*u + theta*kappa.*(T.*xmin-(1/a).*log((1-g.*exp(-T.*d))./(1-g))); 
vector = exp(C + D.*nu0 + 1i.*u.*log(S0)); 
end

