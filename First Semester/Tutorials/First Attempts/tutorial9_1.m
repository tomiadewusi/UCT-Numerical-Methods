clear
clc
rng(0)

T = 1;
S0 = 100;
sigma = 0.4;
r = 0.1;
K = 50;
mu = -2;

d1  = (log(S0/K) + (r + 0.5*sigma^2)*T) /(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);

BSput = exp(-r*T)*K*normcdf(-d2) - S0*normcdf(-d1);

mat1 = zeros(50,3);
mat2 = zeros(50,3);
for idx = 1:50
    n = 1000*idx;
    X = randn(n,1);
    fx = exp(-r*T)*max(K - S0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T).*X),0);
    put_value = mean(fx);
    confidence = 3*std(fx)/sqrt(n);
    lower = BSput - confidence;
    upper = BSput + confidence;
    mat1(idx,:) = [put_value lower upper];
    
    X = X + mu; 
    fx = exp(-r*T)*max(K - S0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T).*X),0);
    wx = normpdf(X,0,1);
    vx = normpdf(X,mu,1);
    sample = fx.*wx./vx;
    put_value = mean(sample);
    confidence = 3*std(sample)/sqrt(n);
    lower = BSput - confidence;
    upper = BSput + confidence;
    
    mat2(idx,:) = [put_value lower upper];
end

figure()
hold on
plot([0 50],[BSput BSput], 'k-')
plot(mat2(:,1), 'r.')
plot(mat2(:,2:3), 'r--')
plot(mat1(:,1), 'b.')
plot(mat1(:,2:3), 'b--')

hold off 
