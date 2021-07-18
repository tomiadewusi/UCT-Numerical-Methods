clear
clc
rng(0)


% Trying replace dec2base
b = 2;
n = 5; 
ix = 1;
cols = ceil(log(n)/log(2));
% this function spits the value in the reverse order as dec2base

mat = zeros(n,cols); 
for jx = 1:n
    a = zeros(1,cols);
    num = jx; 
    while 1
        if idivide(int16(num).',int16(b)) == 0
            a(ix) = mod(num,b);
            
            break
        else
            a(ix) = mod(num,b);
            ix = ix + 1;
            num = idivide(int16(num).',int16(b));
        end
    end
    a
    mat(jx,:) = a; 
end