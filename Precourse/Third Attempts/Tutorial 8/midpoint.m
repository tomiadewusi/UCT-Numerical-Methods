function value = midpoint(f,a,b,h)
	points = a:h:b-h;
	value = h.*(sum(f(points + h/2))); 
end 