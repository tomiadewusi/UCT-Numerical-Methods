clear
clc
rng(0)
% Computing price of the option at inception
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing the intial values
S0 = [35; 65];
sigma1 = 0.25; sigma2 = 0.2;
mu1 = 0.12; mu2 = 0.15;
Sigma = [1 0.6; 0.6 1];
r = 0.1;
T = 1/12;
N = 300;
deltat = T/N;
K = 100;
L = chol(Sigma,'lower');
n = 50001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We now have a Normalised Sobol Sequence of length n
tic
ZS = norminv([Sobol(gendirnums(7,[1 3],16),n);...
              Sobol(gendirnums(11,[1 1 5],16),n)]);
toc
ZS(:,1) = [];
ZS = L*ZS;
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ST = zeros(2,50000);
ST(1,:) = S0(1)*exp((r-0.5*sigma1^2)*T + sigma1*sqrt(T).*ZS(1,:));
ST(2,:) = S0(2)*exp((r-0.5*sigma2^2)*T + sigma2*sqrt(T).*ZS(2,:));

% Price of the option under the Risk Neutral measure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Price = mean(exp(-r*T)*max(sum(ST) - K, 0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num = 1000;
PnL = zeros(num,1);
for jx = 1:num
    % First we need the paths
    Z = L*randn(2,N);
    Si1 = [S0(1) S0(1)*exp(cumsum((mu1-0.5*sigma1.^2)*deltat...
        + sigma1*sqrt(deltat).*Z(1,:),2))];
    Si2 = [S0(2) S0(2)*exp(cumsum((mu2-0.5*sigma2.^2)*deltat...
        + sigma2*sqrt(deltat).*Z(2,:),2))];
    
    % Calculating the hi values
    Zmat1 = repmat(ZS(1,1:5000),N,1); 
    Zmat2 = repmat(ZS(2,1:5000),N,1); 
    time_index = (T - deltat*(0:299))'; 
    Tmat = repmat(time_index,1,5000); 
    
    mat1 = Si1(1:N)'.*exp((r-0.5*sigma1^2).*Tmat+sigma1*sqrt(Tmat).*Zmat1);
    mat2 = Si2(1:N)'.*exp((r-0.5*sigma2^2)*Tmat+sigma2*sqrt(Tmat).*Zmat2);
    indicator_discounted = exp(-r*T).*(mat1 + mat2 > K);
    hi1 = mean(indicator_discounted.*mat1./Si1(1:N)',2);
    hi2 = mean(indicator_discounted.*mat2./Si2(1:N)',2);
    
    
    bi = zeros(N,1);
    bi(1) = Price - hi1(1)*Si1(1) - hi2(1)*Si2(1);
    for ix = 2:N
        bi(ix) = bi(ix-1)*exp(r*deltat) - (hi1(ix)- hi1(ix-1))*Si1(ix)...
            -(hi2(ix)- hi2(ix-1))*Si2(ix);
    end
    PnL(jx) = bi(N)*exp(r*deltat) + hi1(N)*Si1(N+1) + hi2(N)*Si2(N+1) -...
        max(Si1(N+1) + Si2(N+1) - K, 0); 
end
toc
figure()
histogram(PnL,50)
% xlim([-0.8 0.8]);
title('Hedging PnL')

