clear
clc
rng(0)

n = 10000;
X = zeros(n,1);
count = 0;
ix = 1;

while ix <= n
    count = count + 1;
    % Generate a pair of indepedent random variables
    U = rand(1,2);
    u1 = U(1);
    u2 = U(2);
    Y = -log(u1);
    if u2 <= exp(-0.5*(Y-1)^2)
        V = rand(1,1);
        if V < 0.5
            X(ix) = -Y;
            ix = ix + 1;
        else
            X(ix) = Y;
            ix = ix + 1;
        end
    end
end

wastage = n/count; 
X(1:5)

figure() 
subplot(1,2,1) 
hist(X) 
subplot(1,2,2) 
hist(randn(10000,1) ) 

figure()
x = linspace(0,6,1000); 
fx = (2/sqrt(2*pi))*exp(-0.5.*x.^2); 
gx = exp(-x); 
c = sqrt(2*exp(1)/pi); 
subplot(1,2,1) 
hold on 
plot(x,fx,'k-') 
plot(x,c.*gx,'b--')
label1 = sprintf('Density of the absolute value \n of a standard normal number'); 
label2 = 'Scaled exponential density'; 
legend(label1, label2)
hold off 
subplot(1,2,2) 
plot(x, fx./(c.*gx), 'b-') 

