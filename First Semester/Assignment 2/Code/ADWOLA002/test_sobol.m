clear
clc
rng(0)

n = 20001;

tic
% Sobol Sequence of length n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ZS = norminv([Sobol(gendirnums(7,[1 3],16),n);...
    Sobol(gendirnums(11,[1 1 5],16),n)]);
ZS1 = norminv([Sobol(gendirnums(7,[1 3],16),n);...
    Sobol(gendirnums(11,[1 1 5],16),n)]);
toc

tic
temp = randn(2,100000,2);
Z1 = temp(1,:,1);
Z2 = temp(2,:,1);
Z3 = temp(1,:,2);
Z4 = temp(2,:,2);
toc
figure() 
scatter(Z4,Z3,'.') 