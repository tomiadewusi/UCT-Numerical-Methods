clear 
clc 
%always remember that most of the defaults in matlab are column vectors 

sigma = 1/10; 
mus = 1/10; 
mu1 = 1/5; 
mu2 = 1/5; 
k2 = 5/11; 

vs_handle = @(t,v) virus_strain(sigma, mus, mu1, mu2, k2, t, v); 
vs_handle1 = @(t,v) virus_strain1(sigma, mus, mu1, mu2, k2, t, v); 


times = zeros(1000,1); 
for idx = 1: 1000
    tic
    [t, V] = ode45(vs_handle, [0 700], [1 0.0001 0.0001]'); 
    times(idx) = toc; 
end 

times1 = zeros(1000,1); 
for idx = 1: 1000
    tic 
    [t1, V1] = ode45(vs_handle1, [0 700], [1 0.0001 0.0001]'); 
    times1(idx) = toc;  
end 


% figure()
% plot(repmat(t, 1, 3), V )

figure()
hold on 
plot(cumsum(times),'k-')
plot(cumsum(times1),'b-')
%the blue line represents not using unnessary elementwise multiplication 
hold off 

%the verdict is if there you are dealing with scalars then don't use 
%unnessary elementwise operations 
