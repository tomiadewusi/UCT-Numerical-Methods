function vector = forw(f,a,b,h)
	points = a:h:b;
	vector = (f(points + h) - f(points))./h; 
end 