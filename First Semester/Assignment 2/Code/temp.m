clear
clc
%
% gendirnums(13,[1 3 5],10)


 x = [Sobol(gendirnums(13,[1 3 5],10),1000); 
    Sobol(gendirnums(19,[1 1 3 13],10),1000)] 

% tic 

% for ix = 1:500
%     x = [Sobol(gendirnums(13,[1 3 5],10),1000);...
%         Sobol(gendirnums(19,[1 1 3 13],10),1000)];
% end
% toc 