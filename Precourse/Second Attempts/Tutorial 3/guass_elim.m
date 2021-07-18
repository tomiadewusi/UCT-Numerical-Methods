function value = guass_elim(A,b)
    [m,n] = size(A); 
    [~,nn] = size(b); 
    if m~=n
        error('Matrix A not a square')
    end
    if nn  ~= 1 || mm~=m 
        error('Incorrect dimensions for b') 
    end
    M = [A b]; 
    [rows, ~] = size(M); 
    for idxi = 1:rows 
        M(idxi,:) = M(idxi, :)./M(idxi,idxi); 
        for idxj = 1:rows 
            if idxj ~= idxi 
                M(idxj,:) = M(idxj,:) - M(idxi,:).*M(idxj,idxi); 
            end 
        end 
    end 
   
    value = M(:,end); 
    
            
    

end 