clear 
clc
%rng(1,'v4')

U = rand(10000,1); 
u1 = U(1:5000); 
u2 = U(5001:end); 

z1 = sqrt(-2.*log(u1)).*cos(2.*pi.*u2); 
z2 = sqrt(-2.*log(u1)).*sin(2.*pi.*u2);

N1 = [z1;z2]; 

N2 = randn(10000,1); 

figure(1)
subplot(1,2,1)
hist(N1, -5:0.5:5) 
subplot(1,2,2)
hist(N2, -5:0.5:5)

matrix = corr([z1 z2])