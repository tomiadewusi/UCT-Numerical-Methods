function vector = LCG(a,c,m,x1,n) 
x = zeros(n,1); 
x(1) = x1; 
for idx = 2:n
    x(idx) = mod(a*x(idx-1) + c,m); 
end 
vector = x./m; 

end 