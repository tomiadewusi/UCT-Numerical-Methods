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

Paths = S0.*exp(cumsum((r - 0.5*sigma.^2)*h + sigma.*sqrt(h).*Z )); 
AntiPaths = S0.*exp(cumsum((r - 0.5*sigma.^2)*h + sigma.*sqrt(h).*(-Z)));
FinalPaths = [S0.*ones(1,2*n);  [Paths AntiPaths] ]; 

V = @(S) max(K-S,0) ; 
FinalPayoff = V(FinalPaths(end,:));  
Vold = FinalPayoff; 

phi0 = @(x) ones(size(x)); 
phi1 = @(x) 1 - x; 
phi2 = @(x) 1 - 2.*x + 0.5.*x.^2; 

for ix = N:-1:2
    % Remember that there are three different sizes of Vector
    
    Vnew = exp(-r*h).*Vold; % Going backward in time
    EarlyExercise = V(FinalPaths(ix,:)); % Large 
    PostiveEarlyIndicator = EarlyExercise > 0; % Large 
    X = FinalPaths(ix,:); % Large 
    X = X(PostiveEarlyIndicator); % Medium
    Y = (Vnew(PostiveEarlyIndicator)).'; % Medium
    F = [phi0(X); 
         phi1(X); 
         phi2(X)]; 
    BetaHat = (F*F.')\(F*Y); 
    fBetaX = F.'*BetaHat; % Medium
    if ix == N
        figure()
        hold on 
        scatter(X,Y,'.')
        scatter(X,fBetaX,'r.')
        hold off
    end 
    PostiveEarlyExercise = EarlyExercise(PostiveEarlyIndicator);% Med Size
    EarlyOptimalIndicator = (PostiveEarlyExercise > fBetaX.'); % Med Size
    Y(EarlyOptimalIndicator) = PostiveEarlyExercise(EarlyOptimalIndicator);
    % The values that are being replaced in Y are of size Small
    Vnew(PostiveEarlyIndicator) = Y; 
    Vold = Vnew; 
end 


OptionValue = exp(-r*h)*mean(Vold) 
