function vector = guass(A,b)
    [m,n] = size(A); 
    if m~=n 
        error('Matrix A not a square')
    end 

	[mm,nn] = size(b); 

	if mm~=m | nn ~=1
		error('Incorrect dimensions for b')
	end

	M = [A b]; 

	for idx = 1:m
		M(idx,:) = M(idx,:)./M(idx,idx); 
		for idx1 = 1:m
			if idx1 ~= idx
				M(idx1,:) = M(idx1,:) - M(idx,:).*M(idx1,idx);
			end 
		end
	end 
vector = M(:,end); 
end 
 
