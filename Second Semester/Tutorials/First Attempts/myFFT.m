function xhat = myFFT(x)
%myFFT implements the Fast Fourier Transform algorithm as described in the notes
%   Detailed explanation goes here
N = length(x); 
if N == 1
    xhat = x; 
else 
    ehat = myFFT(x(1:2:end)); 
    ohat = myFFT(x(2:2:end)); 
    m = [0 1:((N/2)-1)] ; 
    xhat = [(ehat + ohat.*exp(-1i*(m.*2.*pi./N)))...
        (ehat - ohat.*exp(-1i*(m.*2.*pi./N)))]; 
end

