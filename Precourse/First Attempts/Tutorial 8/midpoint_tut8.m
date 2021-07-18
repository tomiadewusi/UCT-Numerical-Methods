function value = midpoint_tut8(f, a, b, h)
points = a:h:b;
value = h.*sum( f(points + h./2) ) ; 

end

