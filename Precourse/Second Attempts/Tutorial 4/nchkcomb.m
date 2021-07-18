function value = nchkcomb(n,k)

    if k ==1 
        value = [1:n]'; 
    elseif k == n
        value = 1:n; 
    else 
        step1 =  nchkcomb(n-1,k); % has a width of k 
        step2 = nchkcomb(n-1,k-1); % has width of k-1
        % step 3 
        step2(:,k) = n; 
        % step 4 
        value = [step1;step2]; 
    end      
end 