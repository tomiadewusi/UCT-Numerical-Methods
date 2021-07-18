 clear
clc
rng(0)

S0 = 100;
r = 0.085;
sigma = 0.3;
T = 2;
K = 100;
delta_S0 = 5;
S0_1 = S0 + delta_S0;
S0_2 = S0;
 
d1 = (log(S0/K) + T*(r + 0.5*sigma.^2) )/(sigma*sqrt(T));  
d2 = d1 - sigma*sqrt(T); 
act_call_delta = normcdf(d1,0,1); 

%Need to check the equation for the variance
% For the forward difference formula
disp('Time for the Forward difference formula')
tic
mat = zeros(50,3);
for idx = 1:50
    N = idx*1000;
    random_nums = randn(N, 2);
    Z1 = random_nums(:,1);
    Z2 = random_nums(:,2);
    
    f1 = exp(-r*T)*(max(S0_1*exp((r-0.5*sigma.^2)*T +sigma*sqrt(T).*Z1) -K ,0));
    f2 = exp(-r*T)*(max(S0_2*exp((r-0.5*sigma.^2)*T +sigma*sqrt(T).*Z2) -K ,0));
    
    delta_call = (mean(f1-f2))/delta_S0;
    var_call = (1/delta_S0.^2)*(var(f1)+ var(f2))/N; 
    
    upper = act_call_delta + 3*sqrt(var_call); 
    lower = act_call_delta - 3*sqrt(var_call); 
    mat(idx,:) = [delta_call upper lower] ; 
end
toc 
% For the central difference formula using same estimates
disp('Time for the Central difference formula with common values')
tic
mat2 = zeros(50,3);
for idx = 1:50
    S0_1 = S0 + delta_S0;
    S0_2 = S0 - delta_S0;
    N = idx*1000;
    Z = randn(N, 1);
    
    f1 = exp(-r*T)*(max(S0_1*exp((r-0.5*sigma.^2)*T +sigma*sqrt(T).*Z) -K ,0));
    f2 = exp(-r*T)*(max(S0_2*exp((r-0.5*sigma.^2)*T +sigma*sqrt(T).*Z) -K ,0));
    
    delta_call = (mean(f1-f2))/(2*delta_S0);
    cov_call = cov(f1,f2);
    % This line was what was causing my error - This shows the importance 
    % Of using brackets!
    % MATLAB interprets these two lines as the same code. 
    % var_call = delta_S0.^2*(1/4)*(cov_call(1,1)/N+ cov_call(2,2)/N -2*cov_call(2,1)/N); 
    % var_call = (1/4*(delta_S0.^2))*(cov_call(1,1)/N+ cov_call(2,2)/N -2*cov_call(2,1)/N);
    
     var_call = (1/(2*delta_S0.^2))*(cov_call(1,1)+ cov_call(2,2) -2*cov_call(2,1))/N;
    % which is also the same as this 
    % var_call = (1/(N*((2*delta_S0).^2)))*(cov_call(1,1)+ cov_call(2,2) -2*cov_call(2,1));
    
    upper = act_call_delta + 3*sqrt(var_call); 
    lower = act_call_delta - 3*sqrt(var_call); 
    mat2(idx,:) = [delta_call upper lower] ; 
end
toc 

% Pathwise Estimate of Blach Scholes Delta
disp('Time for Pathwise Estimate')
tic
mat3 = zeros(50,3);
for idx = 1:50
    N = idx*1000;
    Z = randn(N,1);
    ST_dash = S0*exp( (r-0.5*sigma.^2)*T + sigma*sqrt(T).*Z);
    indicator = ST_dash > K;
    sample = exp(-r*T)*((ST_dash./S0).*indicator); 
    delta_call = mean(sample);
    confidence = 3*std(sample)/sqrt(N); 
    lower = act_call_delta - confidence;
    upper = act_call_delta + confidence;
    
    mat3(idx,:) = [delta_call lower upper];
end
toc 

% Likelihood Estimate
disp('Time for Likelihood Estimate')
tic
mat4 = zeros(50,3);
for idx = 1:50
    N = 1000*idx;
    Z = randn(N,1);
    ST = S0*exp((r - 0.5*(sigma.^2))*T + sigma*sqrt(T).*Z);
    sample  = exp(-r*T)*(Z.*max(ST-K,0))./(S0*sigma*sqrt(T)); 
    delta_call = mean(sample);
    confidence = 3*std(sample)/sqrt(N); 
    lower = act_call_delta - confidence;
    upper = act_call_delta + confidence;
    mat4(idx,:) = [delta_call lower upper];
end
toc

figure(1)
subplot(2,2,1)
hold on 
plot(mat(:,1),'b.')
plot(mat(:,[2 3]),'b--')
plot([0,50] ,[act_call_delta act_call_delta],'b')
title('Forward difference')
axis tight
hold off

subplot(2,2,2)
hold on 
plot(mat2(:,1),'b.')
plot(mat2(:,[2 3]),'b--')
plot([0,50] ,[act_call_delta act_call_delta],'b')
%axis([0 50 0.6 0.85])
title('Central difference with common values')
axis tight
hold off                       

subplot(2,2,3)
hold on 
plot(mat3(:,1),'b.')
plot(mat3(:,[2 3]),'b--')
plot([0,50] ,[act_call_delta act_call_delta],'b')
%axis([0 50 0.6 0.85])
title('Pathwise Estimate') 
axis tight
hold off                       

subplot(2,2,4)
hold on 
plot(mat4(:,1),'b.')
plot(mat4(:,[2 3]),'b--')
plot([0,50] ,[act_call_delta act_call_delta],'b')
%axis([0 50 0.6 0.85])
title('Likelihood Estimate') 
axis tight
hold off     

