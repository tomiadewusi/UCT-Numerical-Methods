function matrix = nchkcomb(n,k)
	if k == 1 
		matrix = [1:n].'; 
	elseif k == n
		matrix = 1:n;
	else
		step1 = nchkcomb(n-1,k); 
		step2 = nchkcomb(n-1,k-1);
		%Alternative 
		%step2(:,k) = n;
		%step3 = step2; 
		[row,~] = size(step2); 
		temp = n.*ones(row,1);
		step3 = [step2 temp];   
		matrix = [step1; step3]; 
	end 
end  
 
	