clear 
clc

% constructing the tridiagnol matrix 
A = 2.*eye(10,10) - triu(ones(10,10),1) + triu(ones(10,10),2)...
    - tril(ones(10,10),-1) + tril(ones(10,10), -2) ; 


LHS = det(A); 
k = 4; 

RHS = A(k,k).*det(A(1:k-1,1:k-1)) + det(A(1:k-2,1:k-2)); 

Aupper = A(1:5, 1:5); 
Alower = A(end-4:end, end-4:end); 

G = Aupper*Alower; 
% this is also a tridiagonal matrix 

% deleting the first and fifth rows and the second and fourth columns 
G([1 5],:) = []; 
G(:, [2 4]) = []; 

b = ones(10,1); 

% With the diag function, if you pass through a matrix, it returns the
% diagonals as a vector, but if you pass through a vector, it creates a 
% diagonal matrix 

D = diag(diag(G)); 
L = tril(G,-1); 
U = triu(G, 1); 
H = inv(D+L); 
b = ones(3,1); 

xint = zeros(3,1); 
xold = xint; 

for idx = 1:4
    xnew = H*(b - U*xold); 
    % only doing this to get the memo answer
    % the best way would be 
    % xnew = (D+L)\(b - U*xold); 
    xold = xnew; 
end 

check  = G\b; 

[check xnew] 



    