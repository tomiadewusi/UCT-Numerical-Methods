function value = mat2dec(A,b) 
[n,l] = size(A); 
r = b.^(l-1:-1:0); 
R = repmat(r,n,1); 
value = sum(A.*R,2);

% Making more compact using matrix multiplication 
% Replace lines 4 and 5 with 
% value = A*r'; 
end 