function value = FindRuns(vector) 
% RUNS Finds the number of runs in a sequence of numbers 
lenV = length(vector); 
R  = zeros(6,1); 
runs = 1;
for ix = 1:lenV-1
    if vector(ix) >= vector(ix+1)
        if runs > 6
            runs = 6; 
        end 
        R(runs) = R(runs) + 1;
        runs = 1; 
    else
        runs = runs + 1; 
    end 
end
if vector(lenV-1) < vector(lenV)
    if runs > 6
        runs = 6; 
    end 
    R(runs) = R(runs) + 1; 
else
    R(1) = R(1) + 1; 
end 
value = R; 
end 

% The reason why you can just add the runs is because at the end of the
% last run, you went thorough the else statement i.e the run was
% already taken into account. However since you hit the end of the for
% loop, the value for runs couldn't be added to the running total.
% Therefore, all you have to do is check if runs is too big, and then
% add it to the running total 