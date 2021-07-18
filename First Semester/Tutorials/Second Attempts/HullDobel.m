function HullDobel(a,c,m) 
% HULLDOBEL This is a function to see if the 3 Hull-Dobel conditions have
% been met 

% Condition 1
if gcd(c,m) == 1 
    disp("Condition 1 is met") 
end 

% Condition 2 
prime_factors_m = unique(factor(m)); 
flag = 1; 
for ix = 1:length(prime_factors_m) 
    if mod(m,prime_factors_m(ix)) == 0 
        if mod(a-1,prime_factors_m(ix)) ~= 0
            flag = 0; 
        end 
    end 
end 
if flag == 1 
    disp("Condition 2 is met") 
end 

% Condition 3 
if mod(m,4) == 0 
    if mod(a-1,4) == 0
        disp("Condition 3 is met") 
    end 
end 


end 