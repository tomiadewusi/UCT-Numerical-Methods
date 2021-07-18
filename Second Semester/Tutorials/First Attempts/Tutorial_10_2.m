clear
clc
rng(0)
%--------------------------------------------------------------------------
%% Computing Vasicek Bond Prices
%--------------------------------------------------------------------------
r0 = 0.07;
K = 0.15;
alpha = 0.15;
b = 0.09;
sigmav = 0.02;
deltat = 2;
ti = 0:deltat:40;
t = ti; 
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
BetaHat = ones(n,M+1);
F = repmat(FwdInt,n,1);
for ix = 1:M
    BetaHat(:,ix+1) = BetaHat(:,ix).*(1 + deltat.*F(:,ix));
    if ix < N
        mu_j = cumsum((deltat.*F(:,ix+1:end).*sigma.^2)./ ...
            (1+deltat.*F(:,ix+1:end)),2);
        Zi = randn(n,1);
        F(:,ix+1:end) = F(:,ix+1:end).*exp((mu_j-0.5*sigma.^2).*deltat ...
            +sigma*sqrt(deltat).*Zi);
    end
end
sample1 = deltat*max(F-K,0)./BetaHat(:,2:end);
        H1 = mean(sample1,1);
        H1 = H1(2:end); 
        std1 = std(sample1,1);
        std1 = std1(2:end); 
%--------------------------------------------------------------------------
%% The Predictor Corrector LIBOR Algorithim
%--------------------------------------------------------------------------
M = N;
sigma = 0.2;
n = 100000;
BetaHat = ones(n,M+1);
FPC = repmat(FwdInt,n,1);
FPC_ = FPC;
for ix = 1:M
    BetaHat(:,ix+1) = BetaHat(:,ix).*(1 + deltat.*FPC(:,ix)); 
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
end
sample2 = deltat.*max(FPC-K,0)./BetaHat(:,2:end);
        H2 = mean(sample2,1);
        H2 = H2(2:end); 
        std2 = std(sample2,1);
        std2 = std2(2:end); 
%--------------------------------------------------------------------------
%% Getting the Implied Forward LIBOR Rates
%--------------------------------------------------------------------------
HFwd = (H1(1:end-1) - H1(2:end))./(deltat.*H1(2:end));
HPCFwd = (H2(1:end-1) - H2(2:end))./(deltat.*H2(2:end));
%--------------------------------------------------------------------------
%% Analytical Formula for Caplets
%--------------------------------------------------------------------------
d1 = (log(FwdInt(2:end)/K) + 0.5*sigma.^2.*t(2:end-1))./...
    (sigma.*sqrt(t(2:end-1)));
d2 = d1 - sigma.*sqrt(t(2:end-1));
B_call = @(F,T,B) B.*(F.*normcdf(d1) - K.*normcdf(d2));
Caplet = B_call(FwdInt(2:end),t(2:end-1),deltat.*Bonds(3:end));
%--------------------------------------------------------------------------
%% Plotting
%--------------------------------------------------------------------------
 
figure()
subplot(1,2,1)
hold on
plot(Caplet);
errorbar(H1,std1/100,'k.');
hold off

subplot(1,2,2)
hold on
plot(Caplet);
errorbar(H2,std2/100,'k.');
hold off