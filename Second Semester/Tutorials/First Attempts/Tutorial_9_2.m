% ADWOLA002 - Tomi Adewusi
% Tutorial 9 - Problem 2

clear
clc
rng(0)

r0 = 0.10;
alpha = 0.1;
b = 0.1;
sigma = 0.02;
deltat = 0.5;
tN = 5;
tM = 3;
ti = 0:deltat:tN;
K = 98;
sigmaHat = @(t1,t2) sigma.*exp(-alpha.*(t2 - t1));
%--------------------------------------------------------------------------
% Generating the Bond Prices and the Intial Forward Rate Curve
%--------------------------------------------------------------------------
% Vasiecek Model Parameters
A = @(t1,t2) (1./alpha).*(1-exp(-alpha.*(t2-t1)));
C = @(t1,t2) (((sigma.^2)/(2.*alpha.^2))-b) .* ((t2-t1)-A(t1,t2))...
    -((sigma.^2)/(4*alpha)).*A(t1,t2).^2;

BondPrices = exp(-A(0,ti)*r0 + C(0,ti)) ;
% figure()
% plot(ti, BondPrices,'.')
ForwardRates = -log(BondPrices(2:end)./BondPrices(1:end-1))./deltat;
%--------------------------------------------------------------------------
% Discrete HJM Algorithim
%--------------------------------------------------------------------------
mat1 = zeros(50,2);
for im = 1:50
    n = 1000*im;
    Beta = ones(n,1);
    N = length(ti)-1;
    M = N;
    cashflows = zeros(n,M);
    FWD = repmat(ForwardRates,n,1);
    for ix = 1:6
        r = FWD(:,ix); % Vector
        Beta = Beta.*exp(r.*deltat);
        
        if ix < N
            ti_min1 = ti(ix); % scalar
            tj = ti(ix+1:end-1);  % (ix to N-1) row vector
            % (ix to N-1) row vector
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
    FWD = FWD(:,7:end); %Strip out the old short rates 
    BetaNew = exp(cumsum(FWD.*deltat,2)); % Forward Rates at that time
    cashflows = [5 5 5 105]./BetaNew ;
    sample = max(K - sum(cashflows,2),0)./Beta;
    % Dividing by the old Beta to get you to the current time
    mat1(im,:) = [mean(sample) 3*std(sample)/sqrt(n)];
end

JamshidVal = 0.87513;
upper = JamshidVal + mat1(:,2);
lower = JamshidVal - mat1(:,2);

figure()
hold on
plot(mat1(:,1),'b.')
plot(lower,'--')
plot(upper,'--')
plot([0 50] ,[JamshidVal JamshidVal])
hold off

