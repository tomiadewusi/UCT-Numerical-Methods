clear
clc
rng(0)
%--------------------------------------------------------------------------
%% Computing Vasicek Bond Prices
%--------------------------------------------------------------------------
r0 = 0.07;
alpha = 0.15;
b = 0.09;
sigmav = 0.02;
deltat = 2;
ti = 0:deltat:40;
N = length(ti)-1;
A = @(t1,t2) (1./alpha).*(1-exp(-alpha.*(t2-t1)));
C = @(t1,t2) (((sigmav.^2)/(2.*alpha.^2))-b) .* ((t2-t1)-A(t1,t2))...
    -((sigmav.^2)/(4*alpha)).*A(t1,t2).^2;

Bonds = exp(-A(0,ti)*r0 + C(0,ti)) ;
% plot(ti,Bonds)
%--------------------------------------------------------------------------
%% Calculating the Corresponding LIBOR forward rates at Inception
%--------------------------------------------------------------------------
FwdInt = (Bonds(1:end-1) - Bonds(2:end))./(deltat.*Bonds(2:end));
% plot(0:deltat:38,FwdInt)
%--------------------------------------------------------------------------
%% The Simple Discrete LIBOR Algorithim
%--------------------------------------------------------------------------
M = N;
sigma = 0.2;
n = 100000;
BetaHat = ones(n,1);
cashflows = zeros(n,M);
F = repmat(FwdInt,n,1);
for ix = 1:M
    BetaHat = BetaHat.*(1 + deltat.*F(:,ix));
    if ix < N
        mu_j = cumsum((deltat.*F(:,ix+1:end).*sigma.^2)./ ...
            (1+deltat.*F(:,ix+1:end)),2); 
        Zi = randn(n,1); 
        F(:,ix+1:end) = F(:,ix+1:end).*exp((mu_j-0.5*sigma.^2).*deltat ...
            +sigma*sqrt(deltat).*Zi); 
    end
    cashflows(:,ix) = 1./BetaHat; 
end
H = [1 mean(cashflows,1)]; 
%--------------------------------------------------------------------------
%% The Predictor Corrector LIBOR Algorithim
%--------------------------------------------------------------------------
M = N;
sigma = 0.2;
n = 100000;
BetaHat = ones(n,1);
cashflows = zeros(n,M);
FPC = repmat(FwdInt,n,1);
FPC_ = FPC; 
for ix = 1:M
    BetaHat = BetaHat.*(1 + deltat.*FPC(:,ix));
    if ix < N
        mu_j_init = cumsum((deltat.*FPC(:,ix+1:end).*sigma.^2)./ ...
            (1+deltat.*FPC(:,ix+1:end)),2);
        Zi = randn(n,1); 
        FPC_(:,ix+1:end) = FPC(:,ix+1:end).*exp((mu_j_init-0.5*sigma.^2) ...
                         .*deltat + sigma*sqrt(deltat).*Zi); 
        mu_j_term = cumsum((deltat.*FPC_(:,ix+1:end).*sigma.^2)./ ...
            (1+deltat.*FPC(:,ix+1:end)),2);
        FPC(:,ix+1:end) = FPC(:,ix+1:end) ...
            .*exp(0.5.*(mu_j_init+mu_j_term-sigma.^2).*deltat...
            + sigma*sqrt(deltat).*Zi);
    end
    cashflows(:,ix) = 1./BetaHat; 
end
HPC = [1 mean(cashflows)]; 
%--------------------------------------------------------------------------
%% Getting the Implied Forward LIBOR Rates 
HFwd = (H(1:end-1) - H(2:end))./(deltat.*H(2:end)); 
HPCFwd = (HPC(1:end-1) - HPC(2:end))./(deltat.*HPC(2:end));
%--------------------------------------------------------------------------
%% Plotting
%--------------------------------------------------------------------------
% Plotting Bond Prices
figure()
hold on 
plot(ti,Bonds,'b-')
plot(ti,H,'k.')
plot(ti,HPC,'ko')
title('Bond Prices') 
legend('Closed form','Standard','Predictor-Corrector') 
hold off 


figure() 
hold on 
plot(0:deltat:38,FwdInt,'b-')
plot(0:deltat:38,HFwd,'k.') 
plot(0:deltat:38,HPCFwd,'ko') 
title('LIBOR Forward Rates') 
legend('Closed form','Standard','Predictor-Corrector') 
hold off 