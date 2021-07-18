function simple(f,N,startpoint,endpoint)
	s = linspace(startpoint,endpoint,1000);  
	fn =@(n,s,f) min(floor( 2^n.*f(s))./(2^n),n);
	figure(1)
	hold on 
	plot(s,f(s))
	for n = 1:N
		plot(s,fn(n,s,f))
	end
	hold off
end 