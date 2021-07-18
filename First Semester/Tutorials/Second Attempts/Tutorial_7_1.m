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
deltaS0 = 0.1;
S0h = S0 + deltaS0;

mat = zeros(50,2);
delta = zeros(50,1);
path = zeros(50,2);
for ix = 1:50
    n = 1000*ix ;
    Z = randn(n,N ); % You have to use common random variables 
    Si = [S0*ones(n,1) S0*exp(cumsum((r - 0.5*sigma^2)*h...
        + sigma*sqrt(h).*Z,2 ))];
    
    Sih = [S0h.*ones(n,1) S0h*exp(cumsum((r - 0.5*sigma^2)*h...
        + sigma*sqrt(h).*Z,2 ))];
    m_Si = mean(Si,2); % size is nx1
    m_Sih = mean(Sih,2);
    fx = exp(-r*T)*max(m_Si - K, 0);
    fxh = exp(-r*T)*max(m_Sih - K, 0);
    
    delta(ix) = mean(fxh - fx)/deltaS0; 
    
    fxpath = exp(-r*T).*(m_Si > K).*m_Si./S0;
    
    mat(ix,:)  = [mean(fx) 3*std(fx)/sqrt(n)];
    
    path(ix,:) = [mean(fxpath) 3*std(fxpath)/sqrt(n) ];
end

lower = mat(:,1) - mat(:,2);
upper = mat(:,1) + mat(:,2);

lowerpath = path(:,1) - path(:,2);
upperpath = path(:,1) + path(:,2);


figure()
subplot(1,2,1)
hold on
plot(1000:1000:50000, lower, 'b--')
plot(1000:1000:50000, upper, 'b--')
plot(1000:1000:50000, mat(:,1), 'b.')
hold off

subplot(1,2,2)
hold on
plot(1000:1000:50000, delta, 'b.')
plot(1000:1000:50000, path(:,1), 'ro')
plot(1000:1000:50000, lowerpath, 'b--')
plot(1000:1000:50000, upperpath, 'b--')
hold off