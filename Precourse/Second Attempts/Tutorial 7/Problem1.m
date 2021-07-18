clear
clc 
rng(1,'v4')

n = 10000; 

U = rand(n,1); 
u1 = U(1:5000,1);
u2 = U(5001:end,1);

z1 = sqrt(-2.*log(u1)).*cos(2.*pi.*u2); 
z2 = sqrt(-2.*log(u1)).*sin(2.*pi.*u2); 
N1 = [z1;z2]; 
N2 = randn(n,1); 

figure()
subplot(1,2,1)
histogram(N1,-5:0.5:5)
subplot(1,2,2)
histogram(N2,-5:0.5:5)

covariance = cov(N1,N2); 
correlation = corr([z1 z2]); 
%Pairs of columns 
%Concatenate the two vectors 
