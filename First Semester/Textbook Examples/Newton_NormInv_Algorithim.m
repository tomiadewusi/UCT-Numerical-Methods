clear 
clc 
rng(0) 

% Remember that Z should be between in (-1,1) 
z = 0.3; 
goal = norminv(z)
x0 = intial(z); 

episilon = 1; 
tol = 0.0000000001; 
counter = 0; 
xold = x0; 
while episilon > tol && counter < 1000
    xnew = xold + (z - normcdf(xold))*exp(-0.5*xold.^2 + log(sqrt(2*pi))); 
    episilon = xnew - xold; 
    xold = xnew; 
    counter = counter + 1; 
end 
xnew