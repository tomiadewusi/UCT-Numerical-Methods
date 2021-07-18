clear 
clc
rng(0) 

xmax = 4; 
xmin = 2; 
ymax = 3; 
ymin = 1; 

% Closed form Solution 
f = @(x,y) sin(x.*y); 

actual = integral2(f,xmin,xmax,ymin,ymax);


t = 1:50; 
mat = zeros(length(t),4); 
for idx = t
    
    N = idx*1000; 
    rand_numbers = rand(2,N)'; 
    U = rand_numbers(:,1); 
    V = rand_numbers(:,2); 
    
    sample = 4*sin(2 + 2.*U + 4.*V + 4.*U.*V ); 
    estimate = sum(sample)/N; 
    standard = std(sample); 
    confidence = 3*standard/sqrt(N); 
    upper = actual - confidence; 
    lower = actual + confidence; 
    mat(idx,:) = [estimate actual upper lower]; 

end 
figure()
hold on
plot(t, mat(:,2:end) )
plot(t,mat(:,1),'k.')
