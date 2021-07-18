function [root, counter] = bisect( a, b, tol)
%BISECT This function finds the root of a function using the bisection
%method 
%   Detailed explanation goes here
if f(a)*f(b) > 0
    error('The function must have opposite signs at a and b')
end 

    while 1 
        c = 0.5*(a + b);
        if f(c) == 0 | 0.5*(b-a) < tol 
            root = c; 
            break 
        elseif f(a)*f(c) < 0
            b = c;
        elseif f(a)*f(c) > 0 
            a = c;
        end 
    end 
end 

function x1 = f(x)
x1 = x^3 + 5*x^2 + x - 3;
end