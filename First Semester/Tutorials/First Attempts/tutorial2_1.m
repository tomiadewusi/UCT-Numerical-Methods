clear 
clc
rng(0)

sigma = [1  0.9  0.95;
        0.9  1  0.855;
        0.95 0.855 1];   

L = Cholesky(sigma); 
Z = randn(3,10000); 

X = L*Z; 

figure() 
plot3(X(1,:), X(2,:), X(3,:), 'k.')





