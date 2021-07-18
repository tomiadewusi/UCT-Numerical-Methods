clear
clc
%Summary of all functions used

% Analtytical Solution
% value = exact_soln(sigma,K,eta)

% Crude Monte-Carlo - basic version for comparison 
% value = crude_soln(sigma, K, eta, n)

% Control Variate Monte-Carlo - best so far
% value = control_soln(sigma, K, eta, n);

% End strat Monte-Carlo - slow, might need to fix 
% value = end_strat_soln(sigma, K, eta, n); 

% Latin Hypercube sampling - 3x times slower than control variates
% value = lhs_soln(sigma, K, eta, n); 
% Stratification in more than one direction is 
% very inefficient. We can use eigenvectors to find the optimal direction. 

% Antithetic Variate Monte-Carlo - general accuracy not bad, fast
% value = anti_soln(sigma, K, eta, n); 

% Combination of Antithetic Variate Monte-Carlo and Control Variates
% value = control_anti_soln(sigma, K, eta, n); 
% This is the most stable i.e. least number of errors
% This is the fastest method 
% This least amount of error! best so far!
% Set n to 70000 for call options but this needs to be stress tested
% For calls, 50000 is the goal 
% Set n to 20000 for put options  but stress testing  still needs to
% be done

% Error function - Uses 100 samples on a grid
% if plotting = 1 then the surface plot is shown
% This functions shows the total number of errors and total time taken 
% The absolute value of error as described in the assignment
% [no_errors,tot_time, abs_error,benchmark] =...
% error_fun(testfun,eta,n,plotting);

% Need to how end strat and antithetic variables behave together 

% Functions handles
% f = @myfunction

f1 = @end_strat_soln; 



eta = 1;
plotting = 1; 
n = 4000; 
dim = 32; 
[no1,am1, t1,b1]...
    = error_fun(f1,eta,n,plotting,dim); 




