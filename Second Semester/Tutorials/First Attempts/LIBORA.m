clear
clc
rng(0)

% Part a

deltat = 2;
ti = 0:deltat:20;
N = length(ti) - 1;
M = N;
sigma = 0.25;
n  = 50000;

% Observed Market Prices
Bonds = exp( (-7.*ti + (ti+1).*log(ti+1))./(100) );
% Calculating the intital Forward rates
FwdInt = (Bonds(1:end-1) - Bonds(2:end))./(deltat.*Bonds(2:end));
F = repmat(FwdInt,n,1);

H = ones(1,M+1);  
Beta = 1;
for idx = 1:M
    r = log(1+deltat*F(:,idx))/deltat;
    Beta = Beta.*(1 + deltat *F(:,idx));
    if idx<N
        j = idx:N-1;
        mu_j = cumsum((deltat .* F(:,j+1)*sigma^2)./(1+deltat*F(:,j+1)),2);
        Z = randn(n,1);
        F(:,j+1) = F(:,j+1).*exp((mu_j - 0.5*sigma^2 ).*deltat + sigma*sqrt(deltat).*Z);
    end
    H(idx+1) = mean(1./Beta);
end

figure() 
subplot(1,2,1)
hold on 
plot(ti,Bonds)
plot(ti,H) 
hold off 
subplot(1,2,2) 
fwdC = (H(1:end-1) - H(2:end))./(deltat*H(2:end));
hold on 
plot(ti, [nan fwdC] ) 
plot(ti, [nan FwdInt] ) 