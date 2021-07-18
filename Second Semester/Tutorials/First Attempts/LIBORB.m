clear
clc
rng(0)

% Part a

deltat = 0.25;
ti = 0:deltat:5;
N = length(ti) - 1;
M = N;
sigma = 0.25;
sigmas = 0.35; 
S0 = 75; 
n  = 50000;

% Observed Market Prices
Bonds = exp( (-7.*ti + (ti+1).*log(ti+1))./(100) );
% Calculating the intital Forward rates
FwdInt = (Bonds(1:end-1) - Bonds(2:end))./(deltat.*Bonds(2:end));
F = repmat(FwdInt,n,1);
S = zeros(size(F)); 
S(:,1) = S0.*ones(n,1); 
H = ones(1,M);
stdH = zeros(size(H)); 
Beta = 1;
for idx = 1:M
    r = log(1+deltat*F(:,idx))/deltat;
    S(:,idx+1)= S(:,idx).*exp((r-0.5.*sigmas.^2).*deltat+sigmas.*sqrt(deltat).*randn(n,1));   
    Beta = Beta.*(1 + deltat *F(:,idx));
    if idx<N
        j = idx:N-1;
        mu_j = cumsum((deltat .* F(:,j+1)*sigma^2)./(1+deltat*F(:,j+1)),2);
        Z = randn(n,1);
        F(:,j+1) = F(:,j+1).*exp((mu_j - 0.5*sigma^2 ).*deltat + sigma*sqrt(deltat).*Z);
    end
    sample  = max(S(:,idx+1)-S(:,idx),0)./Beta;
    H(idx) = mean(sample); 
    stdH(idx) = std(sample) ; 
end

d1 =@(idx) (0.5.*deltat.*sigmas.^2 - log(Bonds(idx+1)./Bonds(idx)))./(sigmas.*sqrt(deltat)) ; 
d2 = @(idx) d1(idx) - sigmas.*sqrt(deltat); 
P =@(idx) S0.*( normcdf(d1(idx)) - (Bonds(idx+1)/Bonds(idx)).*normcdf(d2(idx))   ); 

figure() 
plot(P(1:4))

% figure()
% subplot(1,2,1)
% hold on
% plot(t(2:end-1),Caplet);
% errorbar(t(2:end-1),H,std1/100,'k.');
% hold off
% subplot(1,2,2)
% hold on
% plot(t(2:end-1),Caplet);
% errorbar(t(2:end-1),H2,std2/100,'k.');
% hold off