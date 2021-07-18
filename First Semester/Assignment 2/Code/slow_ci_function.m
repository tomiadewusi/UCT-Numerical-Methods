% This is me adding some random changes
% Version 2
% c0 = [flip(dec2base(num,2) - dec2base(0,2)) 0 ]; 
% t = find(c0 == 0); 
% value = t(1); 

% Version 1 
% function value = c(num)
% c0 = flip(d2b(num) - '0');
% c0_length = length(c0);
% for ix1 = 1:c0_length
%     if c0(ix1) == 0
%         break
%     end
% end
% if ix1 == c0_length
%     value = c0_length + 1;
% else
%     value = ix1;
% end