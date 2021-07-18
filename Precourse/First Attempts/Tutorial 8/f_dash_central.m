function value = f_dash_central(f,a,b,h)

points = a:h:b;  
value = (f(points+h)-f(points-h)) ./  (2.*h) ;

end 