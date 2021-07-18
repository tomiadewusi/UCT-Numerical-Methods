clear 
clc 
rng(0)

N = 1000; 
deltat = 1/N; 
t = 0:deltat:1; 
S0bar = [70 100 90]'; 
mubar = [0.4 0.1 0.12]'; 
sigmabar = [0.4 0.22 0.25]'; 
sigmamatrix = [1   0.3  0.95; 
              0.3   1   0.55;
              0.95  0.55  1]; 
          
%it needs to be correlated 
L = chol(sigmamatrix,'lower'); 
          
          

 Z = L*randn(3,N); 
 S = zeros(3,N+1); 
 S(:,1) = S0bar; 
 S(1,2:end) = S0bar(1).*cumprod(   exp( (mubar(1) -0.5.*sigmabar(1).^2)...
     .*deltat + sigmabar(1).*sqrt(deltat).*Z(1,:) )   ); 
 
 S(2,2:end) = S0bar(2).*cumprod(   exp( (mubar(2) -0.5.*sigmabar(2).^2)...
     .*deltat + sigmabar(2).*sqrt(deltat).*Z(2,:) )   ); 
 
 S(3,2:end) = S0bar(3).*cumprod(   exp( (mubar(3) -0.5.*sigmabar(3).^2)...
     .*deltat + sigmabar(3).*sqrt(deltat).*Z(3,:) )   ); 
 
 
 %log returns 
 logS = zeros(3,N); 
 logS(1,1:end) = log(S(1,2:end)./S(1,1:end-1)); 
 logS(2,1:end) = log(S(2,2:end)./S(2,1:end-1)); 
 logS(3,1:end) = log(S(3,2:end)./S(3,1:end-1)); 
%  
%  figure() 
%  subplot(1,2,1)
%  plot(t,S')
%  subplot(1,2,2)
%  plot(t(1:end-1),logS')
%  
 c = corrcoef(logS')
          

          