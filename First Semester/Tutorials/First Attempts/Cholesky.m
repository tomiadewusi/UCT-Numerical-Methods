function matrix = Cholesky(sigma)
    [row,col] = size(sigma); 
    
    %Step 1
    if row ~= col 
        error("Matrix must be a square")
    end 
    
    %Step 2 
    L = zeros(row,col); 
    
    %Step 3 
    for idxi = 1:row
        l = sigma(idxi,idxi) - sum(L(idxi,1:idxi-1).^2); 
        if l > 0 
            L(idxi,idxi) = sqrt(l); 
        else
            error('Matrix must a be postive-definite')
        end 
        
        for idxj = idxi + 1:row
            L(idxj,idxi) = (sigma(idxi,idxj) -...
                sum(L(idxi,1:idxi-1).*L(idxj,1:idxi-1)))./L(idxi,idxi); 
        end 
        
    end
    matrix = L; 
end 