clear 
clc 

n = 1500; 
x = magic(n); 
y = x(1,:); 

%Part a 

tic 
tot1 = zeros(size(y)); 
sum = 0 ; 
for idx = 1:length(y)
	sum = sum  + y(idx)^2; 
	tot1(idx) = sum; 
end 
toc 

tic 
tot2 = cumsum(y.^2); 
toc 

%Part b 

tic 
[row, column] = size(x); 
tot3 = 0; 
for idxi = 1:row
	for idxj = 1:column
		if x(idxi,idxj) > 500 && x(idxi,idxj) < 1000
			tot3 = tot3 + 1; 
		end 
	end 
end 
toc 

tic 
tot4 = sum(sum( x > 500 & x < 1000)); 
toc 

%Part c

tic  
arr1 = []; 
for idxi = 1:row
	for idxj = 1:column
		if x(idxi,idxj) > 495 && x(idxi,idxj) < 500
			arr1 = [arr1 ;x(idxi,idxj) ]; 
		end 
	end 
end 
toc 

tic 
arr2 = x( x > 495 & x < 500); 
toc 

