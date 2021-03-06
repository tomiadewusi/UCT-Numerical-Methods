clear 
clc 
rng(0) 

% Implementing the runs function 

test = LCG(2^16+3,0,2^31,2,5000) ; 
x = FindRuns(test) ; 
test1 = LCG(2^16 + 1, 2^8 + 3, 2^32, 10, 5000); 
figure() 
hist(test1) 
X = FindRuns(test1); 
n = 5000; 
mu = n.*[1/6, 5/24, 11/120, 19/720, 29/5040, 1/840]'; 
S =         (1/n).*[4529.4   9044.9  13568  18091 22615    27892; 
                    9044.9  18097.0  27139  36187 45234    55789;
                    13568   27139.0  40721  54281 67852    83685;
                    18091   36187.0  54281  72414 90470   111580;
                    22615   45234.0  67852  90470 113262  139476; 
                    27892   55789.0  83685 111580 139476  172860;
                    ]; 
% Check if the matrix is symmetric                 
all(S == S',"all")

Q = (X - mu)'*S*(X - mu)

teststat = chi2inv(0.99, 6) 