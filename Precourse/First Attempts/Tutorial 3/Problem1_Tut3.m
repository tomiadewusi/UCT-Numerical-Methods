clear 
clc

A = [1 2 -1;
    -2 -6 4;
    -1 -3 3];

b = [1 -2 1].';



[rowA, columnA] = size(A);
[rowb, columnb] = size(b); 

if rowA ~= columnA
    error('Matrix A not a square')
end

if rowA ~= rowb
    error('Incorrect dimensions for b')
end 

%Gaussian elimination 
M = [A b]; 

for idx = 1:rowA
    M(idx,:) = M(idx,:)./M(idx,idx);
    for idx2 = 1:rowA
        if idx2 ~=idx
            M(idx2,:) = M(idx2,:) - M(idx,:).*M(idx2,idx);
        end
    end        
end 

answer = M(:,end);


A1 = [1 2 5 1;
      4 7 6 2;
      5 1 3 5;
      7 2 8 3;];
  
 b1 = [4 6 7 8].';
 
 M1 = [A1 b1];
 
[rowA1, columnA1] = size(A1);
[rowb1, columnb1] = size(b1);
 
for idx = 1:rowA1
    M1(idx,:) = M1(idx,:)./M1(idx,idx);
    for idx2 = 1:rowA1
        if idx2 ~=idx
            M1(idx2,:) = M1(idx2,:) - M1(idx,:).*M1(idx2,idx);
        end
    end        
end 

M1


answer = M1(:,end);
