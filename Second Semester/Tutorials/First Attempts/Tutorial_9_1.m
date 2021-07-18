clear
clc
rng(0)

r0 = 0.07;
alpha = 0.15;
b = 0.09;
sigma = 0.02;
deltat = 0.5;
ti = 0:deltat:10;

sigmaHat = @(t1,t2) sigma.*exp(-alpha.*(t2 - t1));
%--------------------------------------------------------------------------
%% Generating the Bond Prices and the Intial Forward Rate Curve
%--------------------------------------------------------------------------
% Vasicek Model Parameters
A = @(t1,t2) (1./alpha).*(1-exp(-alpha.*(t2-t1)));
C = @(t1,t2) (((sigma.^2)/(2.*alpha.^2))-b) .* ((t2-t1)-A(t1,t2))...
    -((sigma.^2)/(4*alpha)).*A(t1,t2).^2;

BondPrices = exp(-A(0,ti)*r0 + C(0,ti)) ;
% figure()
% plot(ti, BondPrices,'.')
ForwardRates = -log(BondPrices(2:end)./BondPrices(1:end-1))./deltat; 
%--------------------------------------------------------------------------
%% Discrete HJM Algorithim
%--------------------------------------------------------------------------
n = 100000; 
Beta = ones(n,1);
N = length(ti)-1;
M = N;
cashflows = zeros(n,M); 
FWD = repmat(ForwardRates,n,1);
for ix = 1:M
    r = FWD(:,ix); % Vector
    Beta = Beta.*exp(r.*deltat); 
    if ix < N
        ti_min1 = ti(ix); % scalar
        tj = ti(ix+1:end-1);  % (ix to N-1) row vector
        % (ix to N-1) row vector because ti starts at zero 
        sigma_j = sigmaHat(ti_min1.*ones(size(tj)),tj);
        % (ix to N-1) row vector
        mu_j = 0.5.*(      cumsum(sigma_j.*deltat).^2 ...
                       -[0 cumsum(sigma_j(1:end-1).*deltat).^2]   )...
                          ./deltat;
        Zi = randn(n,1); % (n by 1) column vector
        FWD(:,ix+1:end) = FWD(:,ix+1:end) + ...
            mu_j.*deltat + sigma_j.*sqrt(deltat).*Zi; 
    end
    cashflows(:,ix) = 1./Beta; 
end

H = [1 mean(cashflows)]; 
%--------------------------------------------------------------------------
%% Plotting
%--------------------------------------------------------------------------
figure()
hold on 
plot(ti,H,'.')
plot(ti,BondPrices,'-')
title('Bond Prices') 
hold off 

figure() 
hold on 
helper1 = (-log(BondPrices))./([1 0.5:deltat:10]); 
helper2 = (-log(H) )./([1 0.5:deltat:10]); 
plot(1:20,helper1(2:end),'b-')
plot(1:20,helper2(2:end),'k.')  
title('Yields')
hold off 