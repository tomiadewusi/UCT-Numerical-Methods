clear 
clc 
rng(0) 
% Early exercise boundary 
% FD computed the values at an arbitrary time some values will be greater
% than the exercise value find the point at which those values become equal
% to the exercise value 

% The largest lower bound of option values than your early exercise values
sigma = 0.4; 
r = 0.06; 
K = 40; 
l = 0; 
T = 1; 
Smax = 160; 
Smin = l; 
M = 100;
N = 80; 
theta = 0.5; 

deltaS = (Smax - Smin)/N; 
deltaTau = T/M; 

VT = @(S) max(K - S,0) ; % for the intial condition 
V0 = @(S,t) K.*ones(size(S)); % for the boundary conditions 
% You don't discount as you get the value immediately 
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
    (sigma.^2*((Smin/deltaS) + 1)- r).*U0m; 

LastRow = 0.5*deltaTau*(Smax/deltaS - 1)*...
    (sigma.^2*((Smax/deltaS) - 1)+ r).*UNm;

bm = [FirstRow; 
      zeros(N-3,M+1); 
      LastRow] ; 
 
Udiff = zeros(N-1,M+1); 
U0n = Uinit(Smin + (1:N-1).*deltaS); 
Udiff(:,1) = U0n.'; 

omega = 1.55; 
tol = 1e-5; 
I = eye(N-1,N-1); 
xold =  U0n.';  
Svec = (Smin + deltaS.*(1:N-1)).';  
Boundary = zeros(1,M+1); 
for ix = 2:M+1
    err = 10; 
    while err > tol
    A = (theta.*G+(1-theta).*I); 
    y = (((1-theta).*F+theta.*I)*Udiff(:,ix-1)...
        + (1-theta).*bm(:,ix-1) + theta.*bm(:,ix)) ; 
    D = diag(diag(A)); 
    U = triu(A,1); 
    L = tril(A,-1); 
    SOR = (D +omega.*L)\(omega.*y + ((1-omega).*D - omega.*U)*xold);
    EarlyExercise = max(K - Svec ,0);
    Indicator = EarlyExercise < SOR; % Since Early Exercise is less than 
    % the continuation value you hold the option 
    Boundary(ix) = min(Svec(Indicator)); 
    xnew = max(SOR, EarlyExercise); 
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
S = s;
DiffPrice = Udiff(:,end); 
disp('The at the money price is:') 
DiffPrice(S == K)