clear 
clc 
rng(0) 

% Intial Parameters 
%--------------------------------------------------------------------------
r0 = 0.07; 
alpha = 0.15; 
b = 0.09; 
sigma = 0.07; 
n = 50000; 
N = 20; 
R = zeros(n,N+1); 
R(:,1) = r0.*ones(n,1);  
P = ones(N+1,1); 
%--------------------------------------------------------------------------
% Some Useful Inline Functions
%--------------------------------------------------------------------------
d = (4*b*alpha)./(sigma.^2); 
lambda = @(t1,t2,rt1) (4*alpha.*exp(-alpha.*(t2-t1)).*rt1)./...
    (sigma.^2.*(1-exp(-alpha.*(t2-t1)))); 

rt2 =@(t1,t2,rt1) ((sigma.^2).*(1-exp(-alpha.*(t2-t1)))./(4*alpha))...
    .*ncx2rnd(d,lambda(t1,t2,rt1),n,1); % Equation 17
%--------------------------------------------------------------------------
% Closed Form Solution 
%--------------------------------------------------------------------------
gamma = sqrt(alpha.^2 + 2*sigma.^2); 
A = @(t1,t2) (2.*(exp(gamma.*(t2-t1))-1) )./...
    ((gamma+alpha).*(exp(gamma.*(t2-t1))-1)+2*gamma); 

C = @(t1,t2) (2*alpha*b/sigma.^2).*...
    log( (2.*gamma.*exp(0.5.*(alpha+gamma).*(t2 -t1)))./...
    ((gamma+alpha).*(exp(gamma.*(t2-t1))-1)+2*gamma) );  
%--------------------------------------------------------------------------
for ix = 1:20
    t1 = (ix-1).*ones(n,1); 
    t2 = ix.*ones(n,1); 
    R(:,ix+1) = rt2(t1,t2,R(:,ix));
    P(ix+1) = exp(-A(0,t2(1))*r0 + C(0,t2(1))) ; % Analytical Prices
    
end 
%--------------------------------------------------------------------------
AveRate = mean(R,1);  
Y = 0.5.*cumsum(R(:,1:(end-1)) + R(:,2:(end)),2); 
Y1 = cumsum(R(:,1:(end-1)),2);
P_trap = mean(exp(-Y)); 
P_simp = mean(exp(-Y1)); 
figure() 
hold on 
plot(0:20,R(1:20,:),'r-')
plot(0:20, AveRate,'b-')
hold off 
figure() 

hold on 
plot(0:20,P.','b-')
plot(0:20,[1 P_trap],'k.') 
plot(0:20,[1 P_simp],'ko') 
hold off 

figure() 
hold on 
helper1 = (-log(P).')./([1 (1:20)]); 
helper2 = (-log([1 P_trap]) )./([1 (1:20)]); 
helper3 = (-log([1 P_simp]) )./([1 (1:20)]); 
plot(1:20,helper1(2:end),'b-')
plot(1:20,helper2(2:end),'ko') 
plot(1:20,helper3(2:end),'kx') 
hold off 


