clear 
clc 

S0 = 5; 
sigma = 0.4; 
T = 2; 
K = 6; 
r = 0.06; 

L = 10; 
c1 = r*T; 
c2 = (sigma^2)*T; 
c4 = 0; 
a = c1 - L*sqrt(c2 + sqrt(c4)); 
b = c1 + L*sqrt(c2 + sqrt(c4)); 
N = 50; 
n = 0:(N-1); 

% Characteristic Function of Log Stock Price Scaled
% Supports Vectorised Input
charST = @(u) exp(1i.*u.*( log(S0./K)+(r-0.5*sigma.^2).*T)...
    -(0.5*sigma.^2).*T.*u.^2); 

chi_n = @(c,d)(               cos(n.*pi.*((d-a)./(b-a))).*exp(d)...
                             -cos(n.*pi.*((c-a)./(b-a))).*exp(c)...
             +(n.*pi./(b-a)).*sin(n.*pi.*((d-a)./(b-a))).*exp(d)...
             -(n.*pi./(b-a)).*sin(n.*pi.*((c-a)./(b-a))).*exp(c))./...
    (1 + (n.*pi./(b-a)).^2)   ; 

nhelp = n(2:end); 
phi_n = @(c,d) [(d-c),((b-a)./(nhelp.*pi).*sin(nhelp.*pi.*((d-a)./(b-a)))...
                      -(b-a)./(nhelp.*pi).*sin(nhelp.*pi.*((c-a)./(b-a))))]; 


vn = (2./(b-a)).*K.*(chi_n(0,b) - phi_n(0,b)); 
Vsum = real(charST((n.*pi)/(b-a)).*exp(-1i.*n.*pi.*(a/(b-a)))).*vn;
Vsum(1) = Vsum(1).*0.5; 
V = exp(-r*T).*cumsum(Vsum); 

% Black-Scholes Analytical formula 
d1 =@(K) (log(S0./K) + (r + 0.5*sigma.^2).*T)./(sigma.*sqrt(T)); 
d2 =@(K) (log(S0./K) + (r - 0.5*sigma.^2)*T)./(sigma.*sqrt(T));
BS =@(K) (S0.*normcdf(d1(K)) - K.*exp(-r*T).*normcdf(d2(K))); 

BS_price = BS(K); 
Dsum = real(charST((n.*pi)/(b-a)).*...
    ((1i.*n.*pi)./(S0.*(b-a))).*exp(-1i.*n.*pi.*(a/(b-a)))).*vn;
Dsum(1) = Dsum(1).*0.5; 
D = exp(-r*T).*cumsum(Dsum); 


n = n +1; 
figure() 
subplot(1,2,1)
hold on 
plot(n,V)
[V(end) BS_price]
plot([0 50], [BS_price BS_price])
legend('Fourier Price','BS Price')
ylabel('Price') 
xlabel('N') 
hold off 
subplot(1,2,2) 
plot(n,log10(abs(BS_price - V)))
ylabel('Log Absolute Error') 
xlabel('N')

figure() 
subplot(1,2,1)
hold on 
plot(n,D)
plot([0 50], normcdf(d1(K)).*ones(1,2))
legend('Fourier Delta','BS Delta')
ylabel('Delta') 
xlabel('N') 
hold off 
subplot(1,2,2) 
plot((n),log10(abs(normcdf(d1(K)) - D)))
ylabel('Log Absolute Error') 
xlabel('N')

