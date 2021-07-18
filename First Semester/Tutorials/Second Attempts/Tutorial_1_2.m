clear 
clc 

U = LCD(1103515245,12345,2^32,1,10000);
lambda = 0.44;
k = 7; 
X = lambda*((-log(U)).^(1/k)); 
figure() 
hist(X,20) 
mean(X) 
std(X) 