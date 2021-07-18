clear 
clc


f =  @(x) (    -exp(-sqrt(1-x.^2))  +  exp(sqrt(1-x.^2))     ).*(  exp(-0.5*(x.^2))   )./(  sqrt(2*pi)    ) ;
value = integral(f,-1,1); 
value 