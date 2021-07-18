clear
clc
rng(0);
 
r_0 = 0.07;
alpha = 0.15;
b = 0.09;
sigma_v = 0.02;
delta_t = 2;
t = 0:delta_t:40;
K = 0.15;
 
%Vasicek Bond Price functions
A = @(t1,t2) 1/alpha *(1 - exp(-alpha*(t2 - t1)));
C = @(t1,t2) (sigma_v.^2 ./(2.*alpha.^2) - b) .* ( (t2-t1) - A(t1,t2)) - sigma_v.^2 ./(4*alpha) .* A(t1,t2).^2;
B = @(t1,t2) exp(-A(t1,t2).*r_0 + C(t1,t2));
 
%Vasicek Bond Price
BondPrice_vasicek = B(0,t);
 
%LIBOR forward rates
F_vasicek = (B(0,t(1:end-1))-B(0,t(2:end)))./(delta_t*B(0,t(2:end)));
 
sigma = 0.2; %constant forward rate volatility
n = 100000;
N = 20;
M = N;
 
%Closed form
d1 = (log(F_vasicek(2:end)/K) + 0.5*sigma.^2.*t(2:end-1))./(sigma.*sqrt(t(2:end-1)));
d2 = d1 - sigma.*sqrt(t(2:end-1));
B_call = @(F,T,B) B.*(F.*normcdf(d1) - K.*normcdf(d2));
Caplet = B_call(F_vasicek(2:end),t(2:end-1),delta_t.*BondPrice_vasicek(3:end));
 
F = repmat(F_vasicek,n,1);
H = zeros(1,N-1);
std1 = zeros(1,N-1);
std2 = zeros(1,N-1);
 
%Algorithm
F_ = repmat(F_vasicek,n,1);
H2 = zeros(1,N-1);
 
Beta_1 = 1;
Beta_2 = 1;
for idx = 1:M
   r = log(1+delta_t*F(:,idx))/delta_t;
   Beta_1 = Beta_1.*(1 + delta_t *F(:,idx));
   Beta_2 = Beta_2.*(1 + delta_t *F_(:,idx));
   if idx<N
       j = idx:N-1;
       mu_j = cumsum((delta_t .* F(:,j+1)*sigma^2)./(1+delta_t*F(:,j+1)),2);
       Z = randn(n,1);
       F(:,j+1) = F(:,j+1).*exp((mu_j - 0.5*sigma^2 ).*delta_t + sigma*sqrt(delta_t).*Z);
       
        mu_j_init = cumsum((delta_t .* F_(:,j+1)*sigma^2)./(1+delta_t*F_(:,j+1)),2);
       
       
       F_tilde = F_(:,j+1).*exp((mu_j_init - 0.5*sigma^2 ).*delta_t + sigma*sqrt(delta_t).*Z);
       
       
       mu_j_term = cumsum((delta_t .* F_tilde*sigma^2)./(1+delta_t*F_tilde),2);
       
       
       F_(:,j+1) = F_(:,j+1).*exp(0.5*(mu_j_init + mu_j_term - sigma^2 ).*delta_t + sigma*sqrt(delta_t).*Z);
   end
   if idx>1
    % In this formulation Beta is based on the idx+1 values, that's why it
    % is ahead!
    ci_1 = delta_t*max(F(:,idx)-K,0); 
    H(idx-1) = mean(ci_1./Beta_1); %H are the discounted cash flows
    std1(idx-1) = std(ci_1./Beta_1);
    ci_2 = delta_t*max(F_(:,idx)-K,0);
    H2(idx-1) = mean(ci_2./Beta_2); %H are the discounted cash flows
    std2(idx-1) = std(ci_2./Beta_2);
   end
    
end
 
%Libor Rates
forward_H = (H(1:end-1) - H(2:end))./(delta_t*H(2:end));
 
 
figure()
subplot(1,2,1)
hold on
plot(t(2:end-1),Caplet);
errorbar(t(2:end-1),H,std1/100,'k.');
hold off
subplot(1,2,2)
hold on
plot(t(2:end-1),Caplet);
errorbar(t(2:end-1),H2,std2/100,'k.');
hold off