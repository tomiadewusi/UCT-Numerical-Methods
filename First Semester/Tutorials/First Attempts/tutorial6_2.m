clear 
clc 
rng(0)

T = 1; 
r = 0.12; 
sigma = 0.45;
K = 50; 
S0 = 100; 
N = 6; % the number of points on the path
h = T/N; 

% Closed form Solution 

muDash = (r -0.5*sigma.^2)*(T/2);
sigmaDash = sqrt( ( (sigma.^2*T)/(6))*( (T/(T+h)) + 1) );

temp1 = S0*exp(muDash + 0.5*sigmaDash.^2); 
temp2 = normcdf((log(S0/K) + muDash )/(sigmaDash) + sigmaDash); 
temp3 = K*normcdf((log(S0/K) + muDash )/(sigmaDash)); 
c = exp(-r*T)*(temp1*temp2 - temp3);

mat = zeros(50,4);
for idx = 1:50
    n =idx*1000; %the number of simulations
    t = 0:h:1; 
    deltat = h; 
    Z = randn(length(t)-1,n);
    Si = [S0.*ones(1,n); S0*exp(cumsum((r-0.5*sigma.^2)*deltat...
        + sigma*sqrt(deltat).*Z))];

    f_geom = exp(-r*T)*max(prod(Si).^(1/(N+1)) - K,0); 
    cHat = mean(f_geom); 
    lower = c - 3*std(f_geom)/sqrt(n); 
    upper = c + 3*std(f_geom)/sqrt(n);
    mat(idx,:) = [cHat c upper lower];
end 
figure()
hold on
plot(mat(:,2:end) )
plot(mat(:,1),'k.')




