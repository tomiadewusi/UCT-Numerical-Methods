clear 
clc 
rng(0) 

sigma = 0.4; 
r = 0.06; 
K = 40; 
l = 20; 
T = 1; 
Smax = 180; 
Smin = l; 
M = 80;
N = 80; 
theta = 0.75; 

deltaS = (Smax - Smin)/N; 
deltaTau = T/M; 

VT = @(S) max(S - K,0) ; % for the intial condition 
V0 = @(S,t) S-S; % for the boundary conditions 
Vinf = @(S,t) S - K.*exp(-r.*(t)); % for the boundary conditions

Uinit = @(S) VT(S); 
U0 = @(l,t) V0(l,t); 
Uinf = @(S,t) Vinf(S,t) ;

T1 = diag(ones(N-2,1),1) - diag(ones(N-2,1),-1); 
T2 = -2.*eye(N-1,N-1) + abs(T1); 

D1 = diag(Smin/deltaS + (1:N-1)); 
D2 = D1.^2; 

F = (1 - r*deltaTau).*eye(N-1,N-1)...
    + 0.5*r       *deltaTau.*D1*T1...
    + 0.5*sigma.^2*deltaTau.*D2*T2; 

G = 2.*eye(N-1,N-1) - F; 

m = 0:M; 
U0m = U0(Smin.*ones(1,M+1), m.*deltaTau); 
UNm = Uinf(Smax.*ones(1,M+1), m.*deltaTau); 

FirstRow = 0.5*deltaTau*((Smin/deltaS) + 1)*...
    (sigma.^2*((Smin/deltaS) + 1)- r).*U0m; 

LastRow = 0.5*deltaTau*(Smax/deltaS - 1)*...
    (sigma.^2*((Smax/deltaS) - 1)+ r).*UNm;

bm = [FirstRow; 
      zeros(N-3,M+1); 
      LastRow] ; 
 
  
Udiff = zeros(N-1,M+1); 
U0n = Uinit(Smin + (1:N-1).*deltaS); 
Udiff(:,1) = U0n.'; 

omega = 1.4; 
tol = 1e-5; 
I = eye(N-1,N-1); 
xold = U0n.';

for ix = 2:M+1
    err = 10; 
    while err > tol
    A = (theta.*G+(1-theta).*I); 
    y = (((1-theta).*F+theta.*I)*Udiff(:,ix-1)...
        + (1-theta).*bm(:,ix-1) + theta.*bm(:,ix)) ; 
    D = diag(diag(A)); 
    U = triu(A,1); 
    L = tril(A,-1); 
    xnew = (D +omega.*L)\(omega.*y + ((1-omega).*D - omega.*U)*xold);
    err = max(abs(xnew - xold));  
    xold = xnew;
    end 
    Udiff(:,ix) = xnew; 
end 

% You include intial conditions but you exclude boundary conditions 
dim = Smax - Smin -1; 
t = deltaTau.*(0:M); 
s = Smin + deltaS.*(1:N-1); 
[Tt,S] = meshgrid(t,s) ; 
figure()
surf(S,Tt,Udiff)
zlabel('U')
xlabel('S') 
ylabel('t') 
title('Diffusion solution') 


% The analytical solution 
S = s; 
mu = (r - 0.5*sigma.^2)/sigma.^2; 
d1 = (log(S/K))./(sigma.*sqrt(T)) + (1 + mu).*sigma.*sqrt(T); 
d2 = log(l.^2./(S.*K))./(sigma.*sqrt(T)) + (1 + mu).*sigma.*sqrt(T); 
C1 = S.*normcdf(d1) - K.*exp(-r*T).*normcdf(d1 - sigma.*sqrt(T)) - ...
    S.*normcdf(d2).*(l./S).^(2*(mu+1)); 
C2 = K.*exp(-r*T)*normcdf(d2 - sigma.*sqrt(T)).*(l./S).^(2*mu); 
C = C1 + C2; 

DiffPrice = Udiff(:,end); 
figure() 
plot(S,C - DiffPrice.' ) 
disp('The at the money price is:') 
DiffPrice(S == K)