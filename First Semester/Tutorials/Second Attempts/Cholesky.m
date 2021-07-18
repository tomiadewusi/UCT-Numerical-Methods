function value = Cholesky(matrix) 
[r,k] = size(matrix); 
if r~=k
    error("Matrix must be square")
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For ix = 1 
L = zeros(size(matrix)); 
l = matrix(1,1);
if l > 0
    L(1,1) = sqrt(l); 
else 
    error("Matrix must be postive-definite")
end 
 for jx = 1 + 1:k
     L(jx,1) = (1/L(1,1))*(matrix(1,jx) ) ;
 end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ix = 2:(k)
    l = matrix(ix,ix) - sum(L(ix , 1:ix-1).^2); 
    if l > 0
        L(ix,ix) = sqrt(l); 
    else 
        error("Matrix must be postive-definite") 
    end 
    for jx = ix + 1:k
        L(jx,ix) = (1/L(ix,ix))*(matrix(ix,jx) - ...
            sum( L(ix , 1:ix-1).*L(jx , 1:ix-1)) ) ; 
    end 
 value = L; 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 



