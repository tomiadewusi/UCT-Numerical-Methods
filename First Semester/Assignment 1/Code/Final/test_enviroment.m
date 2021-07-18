clear
clc


% Error function - Uses 100 samples on a grid
% if plotting = 1 then the surface plot is shown
% This functions shows the total number of errors and total time taken 
% The absolute value of error as described in the assignment
% [no_errors,tot_time, abs_error,benchmark] =...
% error_fun(testfun,eta,n,plotting);

% Functions handles
% f = @myfunction

f1 = @PriceFadeIn; 

eta = -1;
plotting = 1; 
dim = 100;
[no1,a1, t1,b1] = error_fun(f1,eta,plotting,dim); 

 


