clear 
clc 

diff = 1;  
tol = 1e-7; 
counter = 0; 
sum_old = 0; 

%you want to use AND because if one of the conditions fails, then the 
%the boolean will evaluate to false. if you use OR the condition will only
%fail if both of the conditions fail 
while counter < 1000 && abs(diff) > tol
    points = 0:counter;
    sum_new = sum(    (-1).^(points) ./ (2.*points + 1)   ); 
    diff = sum_new - sum_old; 
    sum_old = sum_new; 
    counter = counter + 1; 
    
end 

format long 
sum_new*4

diff1 = 1;  
tol1 = 1e-7; 
counter1 = 0; 
sum_old1 = 0; 

while counter < 1000 || abs(diff1) > tol1
    points = 0:counter1;
    sum_new1 = sum(   1./((4.*points +1).*(4.*points + 3))   ); 
    diff1 = sum_new1 - sum_old1; 
    sum_old1 = sum_new1; 
    counter1 = counter1 + 1; 
    
end 

sum_new1*8
