clear 
clc 

A = toeplitz([100 1:9]); 
b = ones(10,1); 

Dinv = diag( (1/100).*ones(10,1), 0 ); 
%or Dinv = (1/100).* eye(10,10); 

L = tril(A,-1); 
U = triu(A,1); 

x0 = zeros(10,1); 
xold = x0; 


for idx = 1:5
  xnew = Dinv*(b - (L + U)*xold); 
  xold = xnew; 
end 

format long 
xback = A\b; 
[Ltilder, Utilder] = lu(A); 

y = Ltilder\b; 
xlu = Utilder\y;

[xnew xback xlu] 
