clear
clc
rng(0)

T = 2;
S0 = 100;
sigma = 0.4;
r = 0.1;
K = 100;
d = 50;

d1  = (log(S0/K) + (r + 0.5*sigma^2)*T) /(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
BScall = S0*normcdf(d1) - exp(-r*T)*K*normcdf(d2);

mat1 = zeros(50,3);
mat2 = zeros(50,3);
for idx = 1:50
    n = 1000*idx;
    X = randn(n,1);
    fx = exp(-r*T)*max(S0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T).*X) - K,0);
    put_value = mean(fx);
    confidence = 3*std(fx)/sqrt(n);
    lower = BScall - confidence;
    upper = BScall + confidence;
    mat1(idx,:) = [put_value lower upper];
    
    temp = n/d;
    nums = (((1:50) - 1)/d) + (1/d).*rand(temp,d);
    X = norminv(nums);
    fx = exp(-r*T)*max(S0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T).*X) - K,0);
    put_value = (1/d)*sum(mean(fx));
    put_var = (1/d)*sum(var(fx))/n; 
    confidence = 3*sqrt(put_var);
    lower = BScall - confidence;
    upper = BScall + confidence;
    mat2(idx,:) = [put_value lower upper];
    
end


figure()
hold on
plot([0 50], [BScall BScall],'k-')
plot(mat1(:,1), 'b.')
plot(mat1(:,2:3), 'b--')
plot(mat2(:,1), 'r.')
plot(mat2(:,2:3), 'r--')
hold off