function value = intial(x) 
if x >= 0.5
    s = 1; 
else 
    s = -1; 
end 
value = s*sqrt(abs(-1.6.*log(1.0004 - (1 - 2.*x).^2)));
end 