clear 
clc 
rng(0) 

S0 = 100;
T = 0.5;
r = 0.03;

theta = 0.05;
kappa = 3;
sigma = 0.5;
rho = 0;

K = linspace(90,110,21); 
% Black-Scholes Analytical formula 
d1 = @(sigma,idx) (log(S0./K(idx))+(r+0.5.*sigma.^2)*T)./(sigma.*sqrt(T)); 
d2 = @(sigma,idx) (log(S0./K(idx))+(r-0.5.*sigma.^2)*T)./(sigma.*sqrt(T));
BS = @(sigma,Price,idx) S0*normcdf(d1(sigma,idx))...
                        -K(idx)*exp(-r.*T)*normcdf(d2(sigma,idx))...
    - Price;

RhoPrice1 = zeros(21,1); 
RhoPrice2 = zeros(21,1); 
RhoPrice3 = zeros(21,1); 
for ix = 1:21
    RhoPrice1(ix) = HestPrice(theta,-0.4,sigma,kappa,K(ix)) ; 
    RhoPrice2(ix) = HestPrice(theta,rho,sigma,kappa,K(ix)) ; 
    RhoPrice3(ix) = HestPrice(theta,0.4,sigma,kappa,K(ix)) ; 
end 
ThetaPrice1 = zeros(21,1); 
ThetaPrice2 = zeros(21,1); 
ThetaPrice3 = zeros(21,1); 
for ix = 1:21
    ThetaPrice1(ix) = HestPrice(0.04,rho,sigma,kappa,K(ix)) ; 
    ThetaPrice2(ix) = HestPrice(theta,rho,sigma,kappa,K(ix)) ; 
    ThetaPrice3(ix) = HestPrice(0.06,rho,sigma,kappa,K(ix)) ; 
end 
ThetaVol1 = zeros(21,1);  
ThetaVol2 = zeros(21,1); 
ThetaVol3 = zeros(21,1); 
for ix = 1:21
    fun1 =@(sigma) BS(sigma,ThetaPrice1(ix),ix)- sigma;
    fun2 =@(sigma) BS(sigma,ThetaPrice2(ix),ix)- sigma;
    fun3 =@(sigma) BS(sigma,ThetaPrice3(ix),ix)- sigma;
    ThetaVol1(ix) = fzero(fun1,0.23);
    ThetaVol2(ix) = fzero(fun2,0.23);
    ThetaVol3(ix) = fzero(fun3,0.23);
end

RhoVol1 = zeros(21,1);  
RhoVol2 = zeros(21,1); 
RhoVol3 = zeros(21,1); 

for ix = 1:21
    fun1 =@(sigma) BS(sigma,RhoPrice1(ix),ix)- sigma;
    fun2 =@(sigma) BS(sigma,RhoPrice2(ix),ix)- sigma;
    fun3 =@(sigma) BS(sigma,RhoPrice3(ix),ix)- sigma;
    RhoVol1(ix) = fzero(fun1,0.23);
    RhoVol2(ix) = fzero(fun2,0.23);
    RhoVol3(ix) = fzero(fun3,0.23);
end 

figure() 
subplot(1,2,1)
hold on
plot(K,ThetaVol1)
plot(K,ThetaVol2)
plot(K,ThetaVol3)
hold off 

subplot(1,2,2)
hold on
plot(K,RhoVol1)
plot(K,RhoVol2)
plot(K,RhoVol3)
hold off 