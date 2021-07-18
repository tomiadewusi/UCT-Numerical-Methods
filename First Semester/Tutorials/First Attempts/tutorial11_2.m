clear
clc
rng(0)

% Closed form solution
T = 2;
S0 = 100;
sigma = 0.4;
r = 0.1;
K = 110;
X = 10;

d = (log(S0/K) + (r - 0.5*sigma.^2)*T)/(sigma*sqrt(T));
c = X*exp(-r*T)*normcdf(d);
b = 3;

mat1 = zeros(50,3);
mat2 = zeros(50,3);
Zvc = norminv(vanderCorput(50001,b));
Zvc(1) = [];
for ix = 1:50
    n = 1000*ix;
    Z = randn(n,1);
    ST = S0*exp( (r - 0.5*sigma^2)*T + sigma*sqrt(T).*Z );
    sample = exp(-r*T)*X.*(ST>K);
    option_value = mean(sample);
    confidence = 3*std(sample)/sqrt(n);
    lower = c - confidence;
    upper = c + confidence;
    mat1(ix,:) = [option_value,lower, upper];
    
    Z = Zvc(1:n); 
    ST = S0*exp( (r - 0.5*sigma^2)*T + sigma*sqrt(T).*Z );
    sample = exp(-r*T)*X.*(ST>K);
    option_value = mean(sample);
    Dstar = (log(n)/n)*((b-1)/(4*log(b)));
    V1f = max(diff(sample));
    confidence = V1f*Dstar;
    lower = c - confidence;
    upper = c + confidence;
    mat2(ix,:) = [option_value,lower, upper];
end

figure()
hold on
plot([0 50], [c c],'k-')
plot(mat1(:,1),'b.')
plot(mat1(:,2:3),'b--')
plot(mat2(:,1),'r.')
plot(mat2(:,2:3),'r--')
ylim([3.7 4])
hold off