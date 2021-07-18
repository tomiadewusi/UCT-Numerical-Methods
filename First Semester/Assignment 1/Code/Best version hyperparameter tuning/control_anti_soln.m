function value = control_anti_soln(sigma, K, eta)

% Student no: ADWOLA002
% Student name: Tomi Adewusi
% The variance reduction method is a combination of Antithetics and 
% Control Variates with the control variate being the vanilla option

% Fixed Input values
S0 = 100; L = 70; U = 130; r = 0.1; T = 1; N = 25; deltat = T/(N+1);
 
switch eta
    case 1
        if sigma > 0.19 || K < 85 
            if sigma > 0.24 && K < 91
                n = 100000; 
            else
                n = 62000; 
            end 
        else 
            n = 26000; 
        end 
    case -1 
        if sigma > 0.2 && K > 100
            n = 24000; 
        else 
            n = 9500; 
        end
end

Z0 = randn(n+1000,N+1); 

% Calculation of alpha
Z = Z0(1:1000,:);
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
vanilla = exp(-r*T).*(max(eta.*(Si(:,end ) - K),0));
Si(:,end) = [];
fade = vanilla.*mean(Si > L & Si < U,2);

% Taking antithetics into account 
Z = -Z; 
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
vanilla_m = exp(-r*T).*(max(eta.*(Si(:,end ) - K),0));
Si(:,end) = [];
fade_m = vanilla_m.*mean(Si > 70 & Si < 130,2);

cov_fg = cov(0.5*(fade_m + fade),0.5*(vanilla + vanilla_m));
alpha = -cov_fg(1,2)/cov_fg(2,2);

% Generate standard normal variates and compute Antithetic estimates
Z = Z0(1001:end,:);
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
vanilla = exp(-r*T).*(max(eta.*(Si(:,end ) - K),0));
Si(:,end) = [];
fade = vanilla.*mean(Si > L & Si < U,2);

% taking into account antithetics
Z = -Z; 
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
vanilla_m = exp(-r*T).*(max(eta.*(Si(:,end ) - K),0));
Si(:,end) = [];
fade_m = vanilla_m.*mean(Si > L & Si < U,2);

%Closed form solution
d1 = (log(S0/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
BSoption = eta*(S0*normcdf(eta*d1) - exp(-r*T)*K*normcdf(eta*d2));

% Final Value 
value = mean(0.5*(fade_m + fade) + alpha*(0.5*(vanilla + vanilla_m) - BSoption)); 
end 