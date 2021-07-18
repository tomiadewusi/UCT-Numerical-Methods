clear
clc

diff = 1; 
tol = 1e-7;
counter = 0;  
old = 0;

while (abs(diff) > tol && counter ~= 1000)
    points = 0:counter;
    new =  sum(1./(2.*points + 1).*( (-1).^(points) ) )  ;
    diff = new - old;
    old = new;
    counter = counter + 1; 
    new;
    diff;
end 

clear
clc

diff = 10; 
tol = 1e-7;
counter = 1;  
old = 10;

while (abs(diff) > tol && counter ~= 1000)
    points =0:counter;
    new = sum( (1./(4.*points+1)).*(1./(4.*points + 3)) );
    diff = new - old;
    old = new;
    counter = counter + 1; 
    
end

new*8


    
    
    