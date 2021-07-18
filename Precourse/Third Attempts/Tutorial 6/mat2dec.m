function vector = mat2dec(A,b)
	[n,l] = size(A); 
	r = b.^(l-1:-1:0); 
	vector = A*r.'; 
end 