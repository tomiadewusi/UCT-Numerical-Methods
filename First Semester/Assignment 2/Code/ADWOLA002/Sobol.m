function vector = Sobol(dnums,n)
k = length(dnums);
if n > (2^k + 1)
    error('n too large, either decrease n or provide more direction numbers')
end
sobNums = zeros(1,n);
for ix = 2:n
    sobNums(ix) = bitxor(sobNums(ix-1), dnums(c(ix-2)).*2^(k-c(ix-2)));
end
vector = sobNums./(2^k);
end

function value = c(num)
c0 = flip(d2b(num) - '0', 2);
c0_length = length(c0);
for ix1 = 1:c0_length
    if c0(ix1) == 0
        break
    end
end
if ix1 == c0_length
    value = c0_length + 1;
else
    value = ix1;
end
end 

function value = d2b(d)
[~,e]=log2(max(d));
value=char(rem(floor(d*pow2(1-max(1,e):0)),2)+'0'); 
end 