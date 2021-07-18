clear 
clc 
rng(0)

Sigma = [1     0.9    0.95;
         0.9    1     0.855; 
         0.95  0.855  1      ]; 

L = Cholesky(Sigma);  
Z = randn(3,10000); 
X = L*Z; 
figure() 
scatter3(X(1,:), X(2,:), X(3,:),'k.')

% Vectorised way of checking whether the matrix is symmetric 
% all(all(Sigma == Sigma'))     
% all(Sigma == Sigma',[1 2])
% all(Sigma == Sigma','all')