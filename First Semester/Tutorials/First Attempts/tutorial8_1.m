clear
clc
rng(0)

% Option Parameters
T = 1;
S0 = 100;
sigma = 0.45;
r = 0.12;
K = 50;
N = 6;
h = T/N;

% Calculation of alpha
Z = randn(N,1000);
Si = [S0.*ones(1,1000);S0*exp(cumsum((r-0.5*sigma^2)*h + sigma*sqrt(h).*Z))];
fx = exp(-r*T).*max(mean(Si) - K,0);
gx = exp(-r*T).*max(geomean(Si) - K,0);
cov_fg = cov(fx,gx);
alpha = -cov_fg(1,2)/cov_fg(2,2);

% Closed-Form Solution
mu_dash = (r-0.5*sigma^2)*(T/2);
sigma_dash = sqrt(((sigma^2*T)/(6))*(((T)/(T+h))+1));
temp1 = S0*exp(mu_dash+0.5*sigma_dash^2);
temp2 = normcdf(((log(S0/K)+mu_dash)/sigma_dash)+sigma_dash);
temp3 = K*normcdf((log(S0/K)+mu_dash)/(sigma_dash));
closed_geo = exp(-r*T)*(temp1*temp2-temp3);

mat = zeros(50,3);
mat2 = zeros(50,2);
for idx = 1:50
    n = 1000*idx;
    Z = randn(N,n);
    Si = [S0.*ones(1,n);S0*exp(cumsum((r-0.5*sigma^2)*h + sigma*sqrt(h).*Z))];
    fx = exp(-r*T).*max(mean(Si) - K,0);
    gx = exp(-r*T).*max(geomean(Si) - K,0);
    monte_geo = mean(gx);
    monte_arth = mean(fx); 
    control_arth = monte_arth + alpha*(monte_geo - closed_geo);
    mat(idx,:) = [monte_geo monte_arth control_arth];
    
    confidence_arth = 3*std(fx)/sqrt(n);  
    confidence_control = 3*std(fx + alpha.*gx)/sqrt(n); 

    mat2(idx,:) = [confidence_arth confidence_control]; 
end
best_estimate = control_arth; 

lower_arth = best_estimate - mat2(:,1);
upper_arth = best_estimate + mat2(:,1);

lower_control = best_estimate - mat2(:,2); 
upper_control = best_estimate + mat2(:,2); 

% Plotting
figure()
hold on
plot([0 50], [closed_geo closed_geo],'k-')
plot(mat(:,1),'b*')
plot(mat(:,2),'bx')
plot(mat(:,3),'r.')
plot(lower_arth, 'b-')
plot(upper_arth,'b-')
plot(lower_control, 'r-')
plot(upper_control,'r-')

% legend({'Closed-Form Geometric','Monte-Carlo Geometric'...
%     ,'Monte-Carlo Arithmetic','Control Variate Arthimetic'},...
%     'Location','east')
% ylim([47 50.5])
ylabel('Option Price')
hold off



