clear
clc
% Bond prices for T_0 to T_{N+1}
Bstar=[1.000000000000000 0.931733068514009 0.867336864194859 0.807180507606257 0.751390186399840 ...
    0.699924481765438 0.652631479135652 0.609290924427683 0.569644566898682 0.533417438303370 ...
    0.500332340416184 0.470119346466237 0.442521711052162 0.417299242183633 0.394229917259718 ...
    0.373110314043405 0.353755267713142 0.335997045644106 0.319684243642614 0.304680543400063 ...
    0.290863424934260 0.278122895071373];

% Forward Rates derived from the bond Prices
fstar = [0.070000000000000 0.071163908099692 0.071749550745182 ...
    0.071750938523962 0.071287420411181 0.070454753142712 ...
    0.069338979002964 0.068001160258291 0.066491256944801 ...
    0.064864780213725 0.063149485640679 0.061384231031525 ...
    0.059591220846755 0.057770455086370 0.055955240441108 ...
    0.054162230256338 0.052391424532061 0.050648374383400 ...
    0.048927528695231 0.047239989697800 0.045596859621355];


K = linspace(0.02,0.08,53);
exact = zeros(19,53);
monte = zeros(19,53);
times = zeros(19,53);

exact1 = zeros(19,53);
monte1 = zeros(19,53);
times1 = zeros(19,53);
alpha = 0.084987000000019;
sigma = 0.025512300000002;
r0 = 0.07;
Ti = 1;

for ix = 1:19
    for jx = 1:53
        tic
        monte(ix,jx) = PriceCapletAnti(ix,K(jx));
        times(ix,jx) = toc; 
        exact(ix,jx) = CapletAnalytical(Bstar,fstar,ix,K(jx),r0,alpha,sigma); 
    end
end

for ix = 1:19
    for jx = 1:53
        tic
        monte1(ix,jx) = PriceCapletControlAnti(ix,K(jx));
        times1(ix,jx) = toc; 
        exact1(ix,jx) = CapletAnalytical(Bstar,fstar,ix,K(jx),r0,alpha,sigma); 
    end
end
disp('For the Regular Antithetic')
num_larger_three_percent = sum(sum(100.*abs(1-monte./exact) > 3))
num_larger_two_percent = sum(sum(100.*abs(1-monte./exact) > 2))
num_larger_one_percent = sum(sum(100.*abs(1-monte./exact) > 1))
sum2error1 = sum(sum((monte-exact).^2))
time_taken1 = sum(sum(times) )
disp('---------------------------------------------------------------------')
disp('For the Control Variate corrected Antithetic')
num_larger_three_percent = sum(sum(100.*abs(1-monte1./exact1) > 3))
num_larger_two_percent = sum(sum(100.*abs(1-monte1./exact1) > 2))
num_larger_one_percent = sum(sum(100.*abs(1-monte1./exact1) > 1))
sum2error = sum(sum((monte1-exact1).^2))
time_taken = sum(sum(times1) )
disp('---------------------------------------------------------------------')
(sum2error/sum2error1)*(time_taken/time_taken1)