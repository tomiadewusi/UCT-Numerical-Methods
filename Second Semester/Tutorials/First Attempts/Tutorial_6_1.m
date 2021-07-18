clear 
clc 
rng(0) 

S0  = 50; 
sigma = 0.4; 
T = 1; 
r = 0.06; 
K = 50; 
umax = 30; 
N = 100; 

% Black-Scholes Analytical formula 
d1 = (log(S0/K) + (r + 0.5*sigma.^2)*T)/(sigma*sqrt(T)); 
d2 = (log(S0/K) + (r - 0.5*sigma.^2)*T)/(sigma*sqrt(T));
BS =@(phi) phi*(S0*normcdf(phi*d1) - K*exp(-r*T)*normcdf(phi*d2)); 

% Characteristic function pricing 
charST = @(u) exp(1i.*u.*( log(S0)+(r-0.5*sigma.^2).*T)...
    -(0.5*sigma.^2).*T.*u.^2); 

delta_u = umax/N; 
un = ((1:N) - 0.5).*delta_u; 
k = log(K); 
P1 = 0.5 + (1/pi).*sum(    real(   (  exp(-1i.*un.*k).*charST(un-1i)  )...
    ./(  1i.*un.*charST(-1i)  )   ).*delta_u    ); 
P2 = 0.5 + (1/pi).*sum(real((exp(-1i.*un.*k).*charST(un)...
    ./(1i.*un))).*delta_u);

% Call Price 
format long
C_BS = BS(1)
C_CF = S0*P1 - K*exp(-r*T)*P2

P1_int = real(   (  exp(-1i.*un.*k).*charST(un-1i)  )...
    ./(  1i.*un.*charST(-1i)  )   )  ; 
P2_int = real((exp(-1i.*un.*k).*charST(un)...
    ./(1i.*un))) ; 
% The graph shows that the contribution eventually tapers off to zero 
figure() 
hold on 
plot(un,P1_int) 
plot(un,P2_int)
hold off 

% Put Price
P_BS = BS(-1)
P_CF = -S0*(1-P1) + K*exp(-r*T)*(1-P2)

format short 