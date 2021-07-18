clear 
clc 
rng(0) 

% GBM Parameters
S0 = 50; 
sigma = 0.4; 
T = 1; 
r = 0.06; 

% Algorithim Specific Parameters 
N = 2^10; 
deltav = 0.25; 
alpha = 1.5; 

deltak = (2*pi)/(N*deltav); 
b = 0.5*(N-1)*deltak; 
a = N*deltav; 
n = 0:(N-1); 
vn = (n).*deltav; 
m = 0:(N-1); 
km = -b + m.*deltak; 

% Characteristic for the log-stock price s_T driven by GBM
% This function can handle vectorised input
charST = @(u) exp(1i.*u.*( log(S0)+(r-0.5*sigma.^2).*T)...
    -(0.5*sigma.^2).*T.*u.^2); 

% Fourier Transform of the price of the option with a dappening factor 
% applied  
% This function can handle vectorised input
cHat = @(v) (exp(-r*T).*charST(v -(alpha+1).*1i))./...
    (alpha.^2 + alpha - v.^2 + 1i.*(2*alpha + 1).*v); 

xn = exp(1i.*b.*vn).*cHat(vn).*(deltav/3).*(3 + (-1).^(n+1) -(n==0)); 

xmHat = myFFT(xn);
C_km = ((exp(-alpha.*km))./(pi)).*real(xmHat); 

index = (exp(km) <=100) & (exp(km) >=0.002); 
helper = exp(km); 
Kvec = helper(index) ;
figure()
subplot(1,2,1)
plot(Kvec,C_km(index))
ylabel('European Call Price')
xlabel('Strike')

% Black-Scholes Analytical formula 
d1 =@(K) (log(S0./K) + (r + 0.5*sigma.^2).*T)./(sigma.*sqrt(T)); 
d2 =@(K) (log(S0./K) + (r - 0.5*sigma.^2)*T)./(sigma.*sqrt(T));
BS =@(K) (S0.*normcdf(d1(K)) - K.*exp(-r*T).*normcdf(d2(K))); 

log_abs_err = log10(abs(BS(Kvec) - C_km(index))); 
subplot(1,2,2) 
plot(log10(Kvec),log_abs_err)
ylabel('Log Absolute Error') 
xlabel('Log Strike') 