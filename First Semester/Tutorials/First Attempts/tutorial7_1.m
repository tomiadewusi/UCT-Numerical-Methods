clear
clc
rng(0)

T = 1;
S0 = 100;
delta_S0 = 0.1; 
S0_forward = S0 + delta_S0; 
sigma = 0.45;
r = 0.12;
K = 50;
N = 6;
h = T/N;
mat1 = zeros(50,4);
mat2 = zeros(50,3);

for idx = 1:50
    n = 1000*idx;
    Z = randn(N,n);
   
    Si = [S0.*ones(1,n);S0*exp(cumsum((r-0.5*sigma.^2)*h...
        + sigma*sqrt(h).*Z))];
    Si_forward = [S0_forward.*ones(1,n);S0_forward*...
        exp(cumsum((r-0.5*sigma.^2)*h + sigma*sqrt(h).*Z))];
    
    m_Si = mean(Si); % size is 1xn
    f_arth  =  exp(-r*T)*max(m_Si- K, 0);
    f_arth_forward = exp(-r*T)*max(mean(Si_forward)- K, 0);
    delta_forward = mean(f_arth_forward - f_arth)/delta_S0; 
    
    indicator = m_Si > K; 
    f_arth_pathwise = exp(-r*T).*indicator.*m_Si/S0; 
    
    option_value = mean(f_arth);
     
    delta_pathwise = mean(f_arth_pathwise); 
  
    option_confidence = 3*std(f_arth)/sqrt(n);
    upper = option_value + option_confidence;
    lower = option_value - option_confidence;
    
    delta_p_confidence = 3*std(f_arth_pathwise)/sqrt(n); 
    upper_p = delta_pathwise + delta_p_confidence; 
    lower_p = delta_pathwise - delta_p_confidence;
   
    mat1(idx,:) = [option_value upper lower delta_forward]; 
    mat2(idx,:) = [delta_pathwise upper_p lower_p];
end

figure()
subplot(1,2,1)
hold on
plot(mat1(:,1),'b.')
plot(mat1(:,[2 3]),'b--')
title('Asian arthimetic call option')
hold off
%axis tight

subplot(1,2,2)
hold on
plot(mat1(:,4),'k.')
plot(mat2(:,1),'ro') 
plot(mat2(:,[2 3]),'b--')
hold off
title('Forward vs Pathwise Delta Estimates')
%axis tight 



