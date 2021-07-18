
function A = nchkcomb(n,k)
    %rows = factorial(n)/(factorial(k)*factorial(n-k));
    %temp = zeros(rows, k)
    if k == 1
        A = [1:n].';
       
    elseif k == n
        A = [1:n];
     
    else
        step1 = nchkcomb(n-1,k); 
        step2 = nchkcomb(n-1,k-1);
        %in matlab you can have assign a variable even if it doesnt exist
        %yet.This is and alternate step 3 
        %step2(:,k) = n; %step3
        
        %you don't know before hand how many rows are in step 2 
        [rows,~] = size(step2);
        step2 = [step2 n.*ones(rows,1)];
        A = [step1; step2];
        
    end
   return
end
    
