% This is me adding some random changes
clear 
clc 
d = 5; 
[~,e]=log2(max(d)); % How many digits do we need to represent the numbers?
s=char(rem(floor(d*pow2(1-max(1,e):0)),2)+'0')
dec2bin(d)