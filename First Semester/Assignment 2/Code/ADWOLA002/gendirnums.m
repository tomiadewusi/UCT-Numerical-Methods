function vector = gendirnums(p,m0,k)
d = length(m0);
mi = zeros(1,k);
mi(1:d) = m0;
v2d = 2.^(1:d);
a_i = d2b(p) - '0';
a_i(1) = [];

if length(a_i) > d
    error('too few initial direction numbers for corresponding primitive polynomial')
end 

switch d
    case 2
        for ix = d+1:k
            t = flip(mi(ix - d:ix-1),2);
            t = [(t.*a_i.*v2d) mi(ix-d)];
            mi(ix) = bitxor(t(1),bitxor(t(2),t(3)));
        end
    case 3
        for ix = d+1:k
            t = flip(mi(ix - d:ix-1),2);
            t = [(t.*a_i.*v2d) mi(ix-d)];
            mi(ix) = bitxor(t(1),bitxor(t(2),bitxor(t(3),t(4))));
        end
    case 4
        for ix = d+1:k
            t = flip(mi(ix - d:ix-1),2);
            t = [(t.*a_i.*v2d) mi(ix-d)];
            mi(ix) = bitxor(t(1),bitxor(t(2),bitxor(t(3)...
                ,bitxor(t(4),t(5)))));
        end
end

vector = mi;
end

function value = d2b(d)
[~,e]=log2(max(d));
value=char(rem(floor(d*pow2(1-max(1,e):0)),2)+'0'); 
end 