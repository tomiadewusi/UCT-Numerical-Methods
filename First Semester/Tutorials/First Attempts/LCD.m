function vector = LCD(a,c,m,x1,n)
    vector = zeros(n,1); 
    vector(1) = x1; 
    for idx = 2:n
        vector(idx) = mod(a*vector(idx-1)+c,m); 
    end 
    vector = vector./m; 
end 