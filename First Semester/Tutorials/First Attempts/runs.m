function v = runs(s)
n = length(s);
counter = 1;
v = zeros(6,1); 
for idx = 1:n-1
    if s(idx+1) > s(idx)
        counter = counter + 1; 
    else
        if counter >= 7
            counter = 6;
        end
        if idx == n-1
            break 
        end 
         
        v(counter,1) = v(counter,1) + 1;
        counter = 1;
    end
end

if s(end) > s(end-1) 
    if counter >= 6
        v(6,1) = v(6,1) + 1; 
    else 
        v(counter+1,1) = v(counter+1) + 1; 
    end  
else 
    v(1,1) = v(1,1) + 1; 
    v(counter,1) = v(counter,1) + 1;
end 


end 
