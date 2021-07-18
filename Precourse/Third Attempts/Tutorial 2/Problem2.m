clear 
clc
 
A = 2*eye(10,10) + - triu(ones(10,10),1) + triu(ones(10,10),2)...
-tril(ones(10,10),-1) + tril(ones(10,10),-2); 

values = zeros(10,1);
for k = 3:10
  RHS = A(k,k)*det(A(1:k-1,1:k-1)) + det(A(1:k-2,1:k-2)); 
  values(k) = RHS; 
end 

values([1:2]) = []; 
compare = [repmat(det(A), length(values),1) values];   

upper = A(1:5,1:5); 
lower = A(end-4:end,end-4:end); 

G = upper*lower; 
G(:,[2 4]) = []; 
G([1 5],:) = []; 

b = ones(3,1); 
D = diag(diag(G),0); 
L = tril(G,-1); 
U = triu(G,1); 

%its always good practice to save your intial guess so you can check later
x0 = zeros(3,1);
xold = x0; 
H = inv(D + L); 
b = ones(3,1);

for idx = 1:4
xnew = H*(b - U*xold);
xold = xnew; 
  
end 

format long 
[xnew G\b] 