function value = f_dash_forward(f,a,b,h)

points = a:h:b;  
value = (f(points+h)-f(points)) ./  (h) ; 

end

