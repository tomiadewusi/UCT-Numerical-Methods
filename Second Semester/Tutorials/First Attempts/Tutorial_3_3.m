% ADWOLA002 - Tomi Adewusi 
% American Options Homework

% 1. Find where the continuation value < the early exercise
% value i.e. determine where the option was exercised. 
% 2. Get all the stock values that meet this condition 
% 3. Find the largest stock value that meets this condition 
% Early exercise boundary 
% FD computed the values at an arbitrary time some values will be greater
% than the exercise value find the point at which those values become equal
% to the exercise value 

% The largest lower bound of option values than your early exercise values

clear 
clc 
rng(0) 
sigma = 0.4; 
r = 0.06; 
K = 40; 
T = 1; 
S0 = K; 
N = 50; 
h = T/N; 
n = 50000; 
Z = randn(N,n); 
% Each column is a sample 
% Each row is a timpoint 
Paths = S0.*exp(cumsum((r - 0.5*sigma.^2)*h + sigma.*sqrt(h).*Z )); 
AntiPaths = S0.*exp(cumsum((r - 0.5*sigma.^2)*h + sigma.*sqrt(h).*(-Z)));
FinalPaths = [S0.*ones(1,2*n);  [Paths AntiPaths] ]; 

V = @(S) max(K-S,0) ; 
FinalPayoff = V(FinalPaths(end,:));  
Vold = FinalPayoff; 
phi0 = @(x) ones(size(x)); 
phi1 = @(x) 1 - x; 
phi2 = @(x) 1 - 2.*x + 0.5.*x.^2; 
Boundary = zeros(1,N); 
for ix = N:-1:2
    % Remember that there are three different sizes of Vector
    Vnew = exp(-r*h).*Vold; % Going backward in time
    pathslice = FinalPaths(ix,:); 
    EarlyExercise = V(pathslice); % Large 
    PostiveEarlyIndicator = EarlyExercise > 0; % Large 
    X = pathslice; % Large 
    X = X(PostiveEarlyIndicator); % Medium
    Y = (Vnew(PostiveEarlyIndicator)).'; % Medium
    F = [phi0(X); 
         phi1(X); 
         phi2(X)]; 
    BetaHat = (F*F.')\(F*Y); 
    fBetaX = F.'*BetaHat; % Medium
%     if ix == N
%         figure()
%         hold on 
%         scatter(X,Y,'.')
%         scatter(X,fBetaX,'r.')
%         hold off
%     end 
    PostiveEarlyExercise = EarlyExercise(PostiveEarlyIndicator);% Med Size
    EarlyOptimalIndicator = (PostiveEarlyExercise > fBetaX.'); % Med Size
    ps_trunc = pathslice(PostiveEarlyIndicator); 
    if ix > 7 % There aren't enough samples at the earlier timepoints
        Boundary(ix) = max(ps_trunc(EarlyOptimalIndicator));  
    end 
    Y(EarlyOptimalIndicator) = PostiveEarlyExercise(EarlyOptimalIndicator);
    % The values that are being replaced in Y are of size Small
    Vnew(PostiveEarlyIndicator) = Y; 
    Vold = Vnew; 
end 
Boundary = [Boundary K];            
Boundary(Boundary==0)=nan;  % replace 0 elements with NaN
figure() 
plot(linspace(0,T,51),Boundary)
title('Early Exercise boundary under Simulation') 
ylim([26 40])
xlim([0 T])

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
figure() 
Boundary(1) = K; 
% Backward PDE so you have to reverse it
plot(0:deltaTau:T,Boundary(end:-1:1)) 
title('Early Exercise boundary under Finite Difference') 

