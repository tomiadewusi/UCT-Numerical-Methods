function out = HestPrice(theta,rho,sigma,kappa,K)
%HestPrice Summary of this function goes here
%   Detailed explanation goes here

S0 = 100;
T = 0.5;
r = 0.03;
N = 100;
umax = 30;
nu0 = 0.06; 
% Characteristic function pricing
charST =@(u) HestCharST(u,S0,nu0,theta,kappa,sigma,rho,T,r);

delta_u = umax/N;
un = ((1:N) - 0.5).*delta_u;
k = log(K);
P1 = 0.5 + (1/pi).*sum(    real(   (  exp(-1i.*un.*k).*charST(un-1i)  )...
    ./(  1i.*un.*charST(-1i)  )   ).*delta_u    );
P2 = 0.5 + (1/pi).*sum(real((exp(-1i.*un.*k).*charST(un)...
    ./(1i.*un))).*delta_u);
out = S0*P1 - K*exp(-r*T)*P2;

end

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


