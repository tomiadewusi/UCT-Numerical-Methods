clear
clc
rng(0)

sigma = 0.4;
alpha = 0.9;
r = 0.06;
K = 40;
T = 1;
CEVsigma = @(t,S) sigma.*S.^(alpha - 1) ;

Smin = 0;
Smax = 160;
M = 80;
N = 160;
theta = 0.5;

deltaS = (Smax - Smin)/N;
deltaTau = T/M;

VT = @(S) max(K-S,0) ; % for the intial condition
V0 = @(S,t) exp(-r*t).*K.*ones(size(S)); % for the boundary conditions
Vinf = @(S,t) zeros(size(S)); % for the boundary conditions

Uinit = @(S) VT(S);
U0 = @(l,t) V0(l,t);
Uinf = @(S,t) Vinf(S,t) ;
%-------------------------------------------------------------------------
T1 = diag(ones(N-2,1),1) - diag(ones(N-2,1),-1);
T2 = -2.*eye(N-1,N-1) + abs(T1);
D1 = diag(Smin/deltaS + (1:N-1));
D2 = D1.^2;
I = eye(N-1,N-1);
SigmaM = @(m) diag((CEVsigma(T-m*deltaTau,Smin+(1:N-1).*deltaS)).^2);

m_ = 0:M;
U0m = U0(Smin.*ones(1,M+1), m_.*deltaTau);
UNm = Uinf(Smax.*ones(1,M+1), m_.*deltaTau);

U = zeros(N-1,M+1);
U0n = Uinit(Smin + (1:N-1).*deltaS);
U(:,1) = U0n.';

FirstRow = 0.5*deltaTau*((Smin/deltaS) + 1)*...
    ((CEVsigma(T-m_.*deltaTau,Smin+deltaS)).^2*...
    ((Smin/deltaS) + 1)- r).*U0m;

LastRow = 0.5*deltaTau*(Smax/deltaS - 1)*...
    ((CEVsigma(T-m_.*deltaTau,Smax-deltaS)).^2*...
    ((Smax/deltaS) - 1)+ r).*UNm;

bm = [FirstRow;
    zeros(N-3,M+1);
    LastRow] ;

%for loop begins
for m = 2:M+1
    Fm = (1-r*deltaTau).*I+0.5*r*deltaTau*D1*T1...
        +(0.5.*SigmaM(m))*(deltaTau.*D2*T2);
    
    GmPlus1 = (1+r*deltaTau).*I-0.5*r*deltaTau.*D1*T1...
        -(0.5.*SigmaM(m+1))*(deltaTau.*D2*T2);
    
    A = (theta.*GmPlus1+(1-theta).*I);
    y = (((1-theta).*Fm+theta.*I)*U(:,m-1)...
        + (1-theta).*bm(:,m-1) + theta.*bm(:,m)) ;
    U(:,m) = A\y;
end

% You include intial conditions but you exclude boundary conditions 
dim = Smax - Smin -1; 
t = T - deltaTau.*(0:M); 
s = Smin + deltaS.*(1:N-1); 
[Tt,S] = meshgrid(t,s) ; 
figure()
surf(S,Tt,U)
zlabel('U')
xlabel('S') 
ylabel('t') 
title('Diffusion solution') 
DiffPrice = U(:,end); 
format long 
disp('The at the money price is:') 
DiffPrice(s == K) 
disp('The out the money price is:') 
DiffPrice(s == 2*K) 
format short 

kappa = (2*r)./( (sigma.^2).*(1-alpha).*(exp(2*r*(1-alpha).*T)-1) ); 
x = @(S0) kappa.*S0.^(2*(1-alpha)).*(exp(2*r*(1-alpha).*T)); 
y = kappa*K.^(2.*(1-alpha));
z = 2 + 1/(1-alpha);

PutExact = @(S0) -S0.*ncx2cdf(y,z,x(S0))...
    + K.*exp(-r*T).*(1-ncx2cdf(x(S0),z-2,y));

figure()
hold on 
plot(s,DiffPrice,'k.') 
plot(s,PutExact(s),'b-')
hold off 

figure() 
hold on 
plot(s,PutExact(s)-DiffPrice.' )
hold off 