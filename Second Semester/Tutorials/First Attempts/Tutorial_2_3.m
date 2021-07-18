clear 
clc 
rng(0) 

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
theta = 1; 

deltaS = (Smax - Smin)/N; 
deltaTau = T/M; 

VT = @(S) max(S - K,0) ; % for the intial condition 
V0 = @(S,t) S-S; % for the boundary conditions 
Vinf = @(S,t) S-S; % for the boundary conditions

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
    (sigma.^2*((Smin/deltaS) + 1) - r).*U0m; 

LastRow = 0.5*deltaTau*(Smax/deltaS - 1)*...
    (sigma.^2*((Smax/deltaS) - 1) + r).*UNm;

bm = [FirstRow; 
      zeros(N-3,M+1); 
      LastRow] ; 
  
U = zeros(N-1,M+1); 
U0n = Uinit(Smin + (1:N-1).*deltaS); 
U(:,1) = U0n.'; 
I = eye(N-1,N-1); 

for ix = 2:M+1
    A = (theta.*G+(1-theta).*I); 
    y = (((1-theta).*F+theta.*I)*U(:,ix-1)...
        + (1-theta).*bm(:,ix-1) + theta.*bm(:,ix)) ; 
    U(:,ix) = A\y;  
end 
% You include intial conditions but you exclude boundary conditions 
dim = Smax - Smin -1; 
t = deltaTau.*(0:M); 
s = Smin + deltaS.*(1:N-1); 
[Tt,S] = meshgrid(t,s) ; 
figure()
surf(S,Tt,U)
zlabel('U')
xlabel('S') 
ylabel('t') 
title('Diffusion solution') 

% The analytical solution 
% S = s; 
% n = -100:100; 
% mu = (2*r)/(sigma.^2) + 1; 
% d1 = log(s.
