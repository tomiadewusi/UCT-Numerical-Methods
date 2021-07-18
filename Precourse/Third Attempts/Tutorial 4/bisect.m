function  value = bisect(f,a,b,tol)
	if f(a)*f(b) > 0 
		error('The function must have opposite signs at a and b')
	end 

	c = (a + b)/2; 
	while abs(f(c)) > 0 && (b - a)/2 > tol 
		if f(a)*f(c) < 0 
			b = c; 
		elseif f(a)*f(c) > 0 
			a = c; 
		end 
		c = (a + b)/2;
			
	end 
	value = c; 
end