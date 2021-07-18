clear 
clc 
rng(0)

% In general, you can generate N random numbers in the interval (a,b)
% with the formula r = a + (b-a).*rand(N,1)

t = 1:50; 
mat = zeros(length(t),4); 
for idx = t
    N = idx*1000; 
    U = pi.*rand(N,1); %column vector
    value = (pi/N).*sum(sin(U)); 
    confidence = 3*sqrt(    sum((sin(U) - value).^2)/(N-1)   )/(sqrt(N)); 
    upper = 2 - confidence; 
    lower = 2 + confidence; 
    mat(idx,:) = [value 2 upper lower]; 
end 

figure()
hold on
plot(t, mat(:,2:end) )
plot(t,mat(:,1),'k.')