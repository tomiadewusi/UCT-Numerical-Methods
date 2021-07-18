function value = composite(f,a,b,h) 
    points = a:h:b; 
    value = (h./3).*( f(points(1)) + 4*sum(f(points(2:2:end-1))) +...
        2*sum(f(points(3:2:end-1))) + f(points(end)));    
end 