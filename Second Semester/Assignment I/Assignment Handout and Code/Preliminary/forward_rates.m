clear 
clc 
rng(0) 

N=20;
T_i=0:N+1; % Times T_0 to T_{N+1}
% Bond prices for T_0 to T_{N+1}
Bstar=[1.000000000000000 0.931733068514009 0.867336864194859 0.807180507606257 0.751390186399840 ...
0.699924481765438 0.652631479135652 0.609290924427683 0.569644566898682 0.533417438303370 ...
0.500332340416184 0.470119346466237 0.442521711052162 0.417299242183633 0.394229917259718 ...
0.373110314043405 0.353755267713142 0.335997045644106 0.319684243642614 0.304680543400063 ...
0.290863424934260 0.278122895071373];
TBstar = -log(Bstar); 
endpoint = TBstar(end); 
delta = 0.00000000001; 
r0 = 0.07; 
Ti_trunc = T_i(2:end-1); 

F = griddedInterpolant(T_i, -log(Bstar)); 


helper_fwd = (F(Ti_trunc + delta) - F(Ti_trunc - delta))/(2*delta); 
fwd_rates = [r0, helper_fwd].'; 

figure() 
plot(T_i(1:end-1), fwd_rates) 
 
t1 = zeros(1,22); 
Dates = T_i; 
t2 = Dates; 
F = griddedInterpolant(t2,-log(Bstar)); 
Ti_trunc = Dates(2:end-1); 
fstar = [r0, (F(Ti_trunc + delta) - F(Ti_trunc - delta))/(2*delta)]; 

figure() 
plot(fstar.' - fwd_rates) 