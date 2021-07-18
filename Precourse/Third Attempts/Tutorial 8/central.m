function vector = forward(f,a,b,h)
	points = a:h:b;
	vector = (f(points + h) - f(points - h))./(2.*h); 
end 