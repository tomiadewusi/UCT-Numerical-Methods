clear
clc
rng(0)

% Fixed Parameters
T = 1;
S0 = 100;
sigma = 0.45;
r = 0.12;
K = 50;
N = 6;
h = T/N;

% Closed form solution
mu_bar = (r - 0.5*sigma^2)*T/2 ;
sigma_bar = sqrt( ((sigma^2)*T/6)*(T/(T + h) + 1) );

c1 = S0*exp(mu_bar + 0.5*sigma_bar^2) ;
c2 = normcdf((log(S0/K) + mu_bar)/sigma_bar + sigma_bar ) ;
c3 = K*normcdf((log(S0/K) + mu_bar)/sigma_bar) ;
c = exp(-r*T)*(c1*c2 - c3);

mat = zeros(50,2);
for ix = 1:50
    n = 1000*ix ;
    Si = [S0*ones(n,1) S0*exp(cumsum((r - 0.5*sigma^2)*h + sigma*sqrt(h).*randn(n,N ),2 ))];
    fx = exp(-r*T)*max(geomean(Si,2) - K, 0);
    mat(ix,:) = [mean(fx) 3*std(fx)/sqrt(n)];
end

lower = c - mat(:,2); 
upper = c + mat(:,2); 

figure()
hold on 
plot([0 50000], [c c ], 'k-') 
plot(1000:1000:50000, lower, 'b--') 
plot(1000:1000:50000, upper, 'b--')
plot(1000:1000:50000, mat(:,1), 'b.')
hold off