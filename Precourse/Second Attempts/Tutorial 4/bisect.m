function value = bisect(f,a,b,tol)

    if f(a)*f(b) > 0 
       error("a and b need to produce opposite signs")
    end 

    c = (a + b)/2; 

    while abs(f(c)) > tol && (b-a)/2 > tol 
        if f(a)*f(c) < 0 
            b = c; 
        else 
            a = c; 
        end 
    c = (a + b)/2;

    end 
    value = c;  
end 

