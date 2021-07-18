clear 
clc 
rng(0)

X0 = 50; 
Xmax = 3*X0; 
Xmin = 0; 

mu = 0.06; 
sigma = 0.25; 
T = 1; 

M = 50;
N = 75; 
theta = 0.75; 

deltaS = (Xmax - Xmin)/N; 
deltaTau = T/M; 

VT = @(S) (S==X0)./deltaTau ; % for the intial condition 
V0 = @(S,t) S-S; % for the boundary conditions 
Vinf = @(S,t) S-S; % for the boundary conditions

Uinit = @(S) VT(S); 
U0 = @(l,t) V0(l,t); 
Uinf = @(S,t) Vinf(S,t) ;

T1 = diag(ones(N-2,1),1) - diag(ones(N-2,1),-1); 
T2 = -2.*eye(N-1,N-1) + abs(T1); 

A = @(X) mu.*X - 2.*X.*sigma.^2; 
B = @(X) -0.5.*sigma.^2.*(X.^2); 
C = @(X) (mu - sigma.^2).*ones(size(X)); 


X = Xmin + deltaS.*(1:N-1); 

F = diag(1-C(X)).*deltaTau ...
    - (deltaTau./(2.*deltaS)).*diag(A(X)).*T1 ...
    - (deltaTau/deltaS.^2).*diag(B(X)).*T2 ; 
G = 2.*eye(size(F)) - F; 


m = 0:M; 
U0m = U0(Xmin.*ones(1,M+1), m.*deltaTau); 
UNm = Uinf(Xmax.*ones(1,M+1), m.*deltaTau); 

FirstRow = 0.5*deltaTau*((Xmin/deltaS) + 1)*...
    (sigma.^2*((Xmin/deltaS) + 1)- mu).*U0m; 

LastRow = 0.5*deltaTau*(Xmax/deltaS - 1)*...
    (sigma.^2*((Xmax/deltaS) - 1)+ mu).*UNm;

bm =  zeros(N-1,M+1);  
  
U = zeros(N-1,M+1); 
U0n = Uinit(Xmin + (1:N-1).*deltaS); 
U(:,1) = U0n.'; 
I = eye(N-1,N-1); 

for ix = 2:M+1
    A = (theta.*G+(1-theta).*I); 
    y = (((1-theta).*F+theta.*I)*U(:,ix-1)...
        + (1-theta).*bm(:,ix-1) + theta.*bm(:,ix)) ; 
    U(:,ix) = A\y;  
end 
U; 
% You include intial conditions but you exclude boundary conditions 
dim = Xmax - Xmin -1; 
t = deltaTau.*(0:M); 
s = Xmin + deltaS.*(1:N-1); 
[Tt,S] = meshgrid(t,s) ; 
figure()
surf(S,Tt,U)
zlabel('U')
xlabel('S') 
ylabel('t') 
title('Diffusion solution') 


