clear
clc

n = 1500; 
x = magic(n); 
y = x(1,:); 

% Part A 
% Using for loops 

[~,columns] = size(y); 

tic 
total = 0; 
tot1 = zeros(1,columns); 
for idx = 1:columns 
    total = total + y(1,idx)^2; 
    tot1(1,idx) = total; 
    % note that you can also do tot1 = [tot1 total] 
    % but this is computationally less efficient 
end
toc 

% Using vectorisation
tic 
tot2 = cumsum(y.^2); 
toc 

%Comparision 
check1 =  tot2(1,1000); 
check2 =  tot1(1,1000); 

% Part B 
% Using for loops 

[rows, columns] = size(x); 

tic 
tot3 = 0; 
for idxi = 1:rows
    for idxj = 1:columns 
        if x(idxi,idxj) > 500 && x(idxi,idxj) < 1000
            tot3 = tot3 + 1; 
        end 
    end 
end 
toc 
% Vectorisation 

tic
tot4 = sum(sum( x > 500 & x < 1000)); 
toc 

% Part C 
% Using for loops 
tic 
arr1 = []; 
for idxi = 1:rows
    for idxj = 1:columns 
        if x(idxi,idxj) > 495 && x(idxi,idxj) < 500
            arr1 = [arr1 x(idxi,idxj)] ; 
        end 
    end 
end 
toc 

% Vectorisation 

tic 
arr2 = x(x > 495 & x < 500); 
toc 