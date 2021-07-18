clear 
clc 

A = [1 2 -1; 
	-2 -6 4; 
	-1 -3 3];

b = [1 -2 1].';

vector = guass(A,b); 

A1 = [1 2 5 1; 
	  4 7 6 2;
	  5 1 3 5; 
	  7 2 8 3]; 

b1 = [4 6 7 8]';

vector1 = guass(A1,b1); 