function value = vanderCorput(n,b)
% Step 1
m = 0:n-1;
% Step 2 and Step 3
m = dec2base(m,b) - double('0');
% Step 4 
[~, cols] = size(m); 
temp = -((cols-1):-1:0)-1; 
bi = b.^temp; 
% Step 5
value = m*bi.';
end 