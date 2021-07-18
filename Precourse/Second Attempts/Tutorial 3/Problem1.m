clear 
clc

A = [1 2 -1;
    -2 -6 4;
    -1 -3 3]; 

b = [1 -2 1]';

A1 = [1 2 5 1;
      4 7 6 2;
      5 1 3 5;
      7 2 8 3]; 
  
b1 = [4 6 7 8]';

format long
value = guass_elim(A,b); 
value1 = guass_elim(A1,b1); 

