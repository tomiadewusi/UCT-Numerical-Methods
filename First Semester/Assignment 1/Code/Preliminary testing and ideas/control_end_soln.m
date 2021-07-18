function value = control_end_soln(sigma, K, eta, n)
% Fixed Input values
S0 = 100; L = 70; U = 130; r = 0.1; T = 1; N = 25; t = linspace(0,T,N+2); 

% n0  should be 10% of the total number of samples 

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
W = BrownianBridge(t,W0,W100,Z1); 
Si = S0*exp((r - 0.5*sigma^2)*t(2:end-1) + sigma.*W);
vanilla = exp(-r*T).*(max(eta.*(ST - K),0)); 
fade = vanilla.*mean(Si > L & Si < U,2);
% Preparing for Stratication
fade = reshape(fade, n/d,d); 
vanilla = reshape(vanilla, n/d,d);
% Trying to capture as much of the between strata variance as possible
% Trying to sample across all the strata
% Since d is 100, I a using the first 1000 points for my calculation of
% alpha
fade_strata = reshape(fade(1:10,:),[],1); 
vanilla_strata = reshape(vanilla(1:10,:),[],1); 
% Calculation for alpha 
cov_fg = cov(fade_strata,vanilla_strata);
alpha = -cov_fg(1,2)/cov_fg(2,2);
%Closed form solution
d1 = (log(S0/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
BSoption = eta*(S0*normcdf(eta*d1) - exp(-r*T)*K*normcdf(eta*d2));
% Final Value 
value = (1/d)*sum(mean(fade + alpha.*(vanilla - BSoption)));
end 

function W = BrownianBridge(t,W0,WN,Z)

if isempty(Z) 
    W = [];
else

    N = length(t)-1;
    m = floor((N+1)*0.5);

    Wm = (((t(1,end)-t(1,m+1))*W0+(t(1,m+1)-t(1,1))*WN)/(t(1,end)-t(1,1)))...
       + sqrt((t(1,m+1)-t(1,1))*(t(1,end)-t(1,m+1))/(t(1,end)-t(1,1))).*Z(:,m);

    Wl = BrownianBridge(t(1:m+1),W0,Wm,Z(:,1:m-1));
    Wr = BrownianBridge(t(m+1:end),Wm,WN,Z(:,m+1:end));
    W = [Wl Wm Wr]; 
end 
end