clear 
clc

%building the diagnol matrix 
temp = eye(10, 10).*2;

%building the -1 on off diagnol 
temp2 = tril(ones(10, 10), 1) - tril(ones(10, 10), 0);
temp3 = tril(ones(10, 10), -1) - tril(ones(10, 10), -2);

temp4 = -temp2-temp3;

final = temp4 + temp;

k = 4;
final_det = det(final); 

final_k = final(k,k);
det_k_minus_1 = det( final(1:(k-1), 1:(k-1)) ); 
det_k_minus_2 = det( final(1:(k-2), 1:(k-2)) ); 

rhs = final_det;
lhs = final_k.*det_k_minus_1 + det_k_minus_2; 

%this value is very close to zero
rhs - lhs;

%these are same 
UL5x5 = final(1:5, 1:5);
LR5x5 = final(end-4:end, end-4:end);

G = UL5x5*LR5x5;
%G is still a tridiagnol matrix 

G([1 5],:) = [];
G(:,[2 4]) = [];


b = ones(3,1);

D  = diag(diag(G), 0); 

L = tril(G, -1);
U = triu(G, 1);

H = inv(D+L); 

xint = zeros(3, 1); 

x1 = H*(b - U*xint);
x2 = H*(b - U*x1);
x3 = H*(b - U*x2);
x4 = H*(b - U*x3);
x5 = H*(b - U*x4);

x5

xBD = G\b





