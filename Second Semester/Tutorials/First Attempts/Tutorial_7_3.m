clear 
clc 

S0 = 7; 
sigma = 0.04; 
T = 2; 
K = 6; 
X = 10; 
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


vn = X.*(2./(b-a)).*(phi_n(a,0)); 

u = ((n.*pi)./(b-a)); 
Vsum = real(charST(u).*exp(-1i.*u.*a) ).*vn;
Vsum(1) = Vsum(1).*0.5; 
V = exp(-r*T).*cumsum(Vsum); 

u = ((n.*pi)./(b-a)); 
Dsum = real( charST(u).*(-sigma.*T.*(1i.*u + u.^2) )...
        .*exp(-1i.*u.*a) ).*vn; 
Dsum(1) = Dsum(1).*0.5; 
D = exp(-r*T).*cumsum(Dsum); 



d2 = (log(S0./K) + (r - 0.5*sigma.^2)*T)./(sigma.*sqrt(T));
BS_price = X*exp(-r*T)*normcdf(-d2); 
BS_vega = X*exp(-r*T)*normpdf(-d2)*(sqrt(T) + d2/sigma); 
n = n +1; 
figure() 
subplot(1,2,1)
hold on 
plot(n,V)
plot([0 N], [BS_price BS_price])
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
plot([0 N], [BS_vega BS_vega])
legend('Fourier Vega','BS Vega')
ylabel('Price') 
xlabel('N') 
hold off 
subplot(1,2,2) 
plot(n,log10(abs(BS_vega - D)))
ylabel('Log Absolute Error') 
xlabel('N')
