clear
clc

A = toeplitz([100 1:9]); 

%extracting the diagnol elements from the matrix 

Dinv = eye(10, 10).*(1/100);
%returns the elements on and below the kth diagonal of A. starts at k =0 
Alower = tril(A, -1);
Aupper = triu(A, 1);

b = ones(10, 1);

%initial guess - make sure that the matrix is conformable!
xint = zeros(10, 1);

%iteration 1 
x1 = Dinv*(b - (Alower + Aupper)*xint);

%iteration 2 
x2 = Dinv*(b - (Alower + Aupper)*x1);

%iteration 3 
x3 = Dinv*(b - (Alower + Aupper)*x2);

%iteration 4
x4 = Dinv*(b - (Alower + Aupper)*x3);

%iteration 5 
x5 = Dinv*(b - (Alower + Aupper)*x4);

x5

%using the back division operator 
%this is the same as inv(A)*b
xBD = A\b;

%LU decomposition 

[L,U] = lu(A);


y = L\b; 
xLU =U\y;
xBD

