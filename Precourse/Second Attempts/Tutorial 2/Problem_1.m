clear 
clc

A = toeplitz([100 1:9]); 
b = ones(10,1); 

Dinv = eye(10,10).*(1/100);
% Alernative way 
% Dinv = diag([ones(10,1)]).*(1/100)

% Note that when a function requires a vector as an input
% always concatenate - [] - it first just to be safe

% note that lower and triangular matrices should not have the main diagnols
% included by definition 

L = tril(A, -1); 
U = triu(A, 1); 

xint = zeros(10,1);
xold = xint; 

for idx = 1:5
    xnew = Dinv*(b - (L+U)*xold); 
    xold = xnew; 
end

[Ltilder, Utilder] = lu(A);
y = Ltilder\b; 
xlu = Utilder\y; 

format long
[xnew xlu] 









