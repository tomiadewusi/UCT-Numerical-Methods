function value = end_strat_soln(sigma, K, eta, n)
% Fixed Input values
S0 = 100; L = 70; U = 130; r = 0.1; T = 1; N = 25; t = linspace(0,T,N+2); 

d = 100; 
index = 1:d;
Unums = ((index - 1)/d) + (1/d).*rand(n/d,d);
Z = norminv(Unums); % Used for the Wiener process terminal values
% The Brownian Bridge
Z = reshape(Z,[],1);
Z1 = randn(n,N); % generating all at the same time
ST = S0*exp((r - 0.5*sigma^2)*T + sigma*sqrt(T).*Z);
W0 = zeros(n,1);
W100 = ((log(ST./S0) - (r - 0.5*sigma^2)*T)/sigma);
W = BrownianBridgeV(t,W0,W100,Z1); 
Si = S0*exp((r - 0.5*sigma^2)*t(2:end-1) + sigma.*W);
vanilla = exp(-r*T).*(max(eta.*(ST - K),0)); 
fade = vanilla.*mean(Si > L & Si < U,2);
fade = reshape(fade, n/d,d); % so that suitable for Stratication
value = (1/d)*sum(mean(fade));

end 

function W = BrownianBridgeV(t,W0,WN,Z)

if isempty(Z) 
    W = [];
else

    N = length(t)-1;
    m = floor((N+1)*0.5);

    Wm = (((t(1,end)-t(1,m+1))*W0+(t(1,m+1)-t(1,1))*WN)/(t(1,end)-t(1,1)))...
       + sqrt((t(1,m+1)-t(1,1))*(t(1,end)-t(1,m+1))/(t(1,end)-t(1,1))).*Z(:,m);

    Wl = BrownianBridgeV(t(1:m+1),W0,Wm,Z(:,1:m-1));
    Wr = BrownianBridgeV(t(m+1:end),Wm,WN,Z(:,m+1:end));
    W = [Wl Wm Wr]; 
end 
end