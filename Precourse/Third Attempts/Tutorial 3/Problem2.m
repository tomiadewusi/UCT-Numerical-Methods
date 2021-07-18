clear 
clc 

tol = 1e-7; 
diff = 1; 
oldsum = 0; 
counter = 0; 
newsum = 0; 

while counter < 1000 && abs(diff) > tol
	newsum = newsum + (-1)^(counter)*(1/(2*counter + 1)); 
	diff = newsum - oldsum; 
	oldsum = newsum; 
	counter = counter + 1;  
end 

diff1 = 1; 
oldsum1 = 0; 
counter1 = 0; 
newsum1 = 0; 

while counter1 < 1000 && abs(diff1) > tol
	newsum1 = newsum1 + (1/(4*counter1 + 1))*(1/(4*counter1 + 3)); 
	diff1 = newsum1 - oldsum1; 
	oldsum1 = newsum1; 
	counter1 = counter1 + 1;  
end 