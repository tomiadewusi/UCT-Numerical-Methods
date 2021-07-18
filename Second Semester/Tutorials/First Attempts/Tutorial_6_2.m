clear 
clc 
rng(0) 

x = [zeros(1,16) ones(1,32) zeros(1,16)]; 
xhat = myFFT(x); 
figure() 
plot(x)
title('Orginal Sequence') 
figure() 
plot(real(xhat)) 
title('Real Part of x')
figure()
plot(imag(xhat))
title('Imaginary Part of x') 

