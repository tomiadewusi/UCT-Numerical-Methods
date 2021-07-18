clear
clc
rng(0)

% Initial Parameters
r0 = 0.07;
alpha = 0.1;
b = 0.09;
sigma = 0.025;
K = 0.8;
T = 1;
S = 2;

% Useful Inline Functions
% Need to make sure that they can handle vectorised output
%--------------------------------------------------------------------------
A = @(t1,t2) (1./alpha).*(1-exp(-alpha.*(t2-t1)));
C = @(t1,t2) (((sigma.^2)/(2.*alpha.^2))-b) .* ((t2-t1)-A(t1,t2))...
    -((sigma.^2)/(4*alpha)).*A(t1,t2).^2;

mu_r = @(t1,t2,rt1) rt1 + alpha.*(b - rt1).*A(t1,t2);
sigma_r =@(t1,t2) sqrt((0.5.*sigma.^2).*A(2.*t1, 2.*t2));

mu_Y = @(t1,t2,Yt1,rt1) Yt1 + (t2-t1).*b + (rt1-b).*A(t1,t2);
sigma_Y = @(t1,t2) sqrt( (sigma.^2/alpha.^2) .*...
    ((t2 - t1) - A(t1,t2) - 0.5.*alpha.*A(t1,t2).^2 ));

sigma_rY = @(t1,t2) (0.5.*sigma.^2).*A(t1,t2).^2;

rho_rY = @(t1,t2) (sigma_rY(t1,t2))./( sigma_r(t1,t2).*sigma_Y(t1,t2) );
%--------------------------------------------------------------------------
mat1 = zeros(50,2);
for ix = 1:50
    %start the for loop
    n = ix*1000;
    t0 = zeros(n,1);
    t1 = T.*ones(n,1);
    t2 = S.*ones(n,1);
    
    rt0 = r0.*ones(n,1);
    Yt0 = zeros(n,1);
    Z = randn(n,2);
    Z1 = Z(:,1);
    Z2 = Z(:,2);
    % These are the joint realisations for the short rate and corresponding
    % discount factors for time T
    %----------------------------------------------------------------------
    rt1 = mu_r(t0,t1,rt0) + sigma_r(t0,t1).*Z1;
    Yt1 = mu_Y(t0,t1,Yt0,rt0) + sigma_Y(t0,t1).*...
        (rho_rY(t0,t1).*Z1 + sqrt(1-rho_rY(t0,t1).^2).*Z2 );
    BetaT = exp(Yt1);
    
    Z = randn(n,2);
    Z1 = Z(:,1);
    Z2 = Z(:,2);
    
    rt2 = mu_r(t1,t2,rt1) + sigma_r(t1,t2).*Z1;
    Yt2 = mu_Y(t1,t2,Yt1,rt1) + sigma_Y(t1,t2).*...
        (rho_rY(t1,t2).*Z1 + sqrt(1-rho_rY(t1,t2).^2).*Z2 );
    Beta1 = exp(Yt2);
    %----------------------------------------------------------------------
    % Computation of the Bond Price using the bond pricing formula
    BTS = exp(-A(t1,t2).*rt1 + C(t1,t2));
    % From simulation
    BTS1 = exp(-0.5.*(rt1 + rt2));
    %----------------------------------------------------------------------
    % Bond Option Value
    sample = max(BTS-K,0)./BetaT;
    mat1(ix,1) = mean(sample);
    mat1(ix,2) = 3*std(sample)/sqrt(n);
end
%--------------------------------------------------------------------------
% Closed form solution
%--------------------------------------------------------------------------
B0T = exp(-A(0,T).*r0 + C(0,T));
B0S = exp(-A(0,S).*r0 + C(0,S));
sigmaBar = sigma*sqrt( (1-exp(-2*alpha*T))/(2*alpha) ).*A(T,S);
d1 = ( log(B0S/(B0T*K)) + (0.5*sigmaBar.^2))/sigmaBar ;
d2 = d1 - sigmaBar;
B0_closed = B0S*normcdf(d1) - K*B0T*normcdf(d2);

lower = B0_closed - mat1(:,2);
upper = B0_closed + mat1(:,2);
%--------------------------------------------------------------------------
figure()
hold on
plot([0 50], B0_closed.*ones(1,2),'k-')
plot((1:50).', mat1(:,1), 'b.')
plot(lower, 'b--')
plot(upper, 'b--')
hold off
