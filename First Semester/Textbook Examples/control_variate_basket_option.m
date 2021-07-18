clear
clc
rng(0)

sigma_mat = [1, 0.8,0.9;
            0.8, 1, 0.75;
            0.9 0.75, 1];
L = chol(sigma_mat, 'lower');

S0a = 25;
sigma_a = 0.4;
S0b = 55;
sigma_b = 0.25;
S0c = 120;
sigma_c = 0.3;

T = 1.5;
r = 0.1;

% Closed form solution 
d1 = ((r + 0.5*sigma_c^2)*sqrt(T))/(sigma_c);
d2 = d1 - sigma_c*sqrt(T);
% An at the money European Call => K = ST
c0_c = S0c*(normcdf(d1) - exp(-r*T)*normcdf(d2));

mat = zeros(50,3);
for idx = 1:50
    n = 1000*idx;
    ratio_use = 0.1;
    n0 = n*ratio_use;
    
    Z_temp = randn(3,n);
    Z = L*Z_temp;
    
    ST = zeros(3,n);
    ST(1,:) = S0a*exp( (r-0.5*sigma_a^2)*T + sigma_a*sqrt(T).*Z(1,:) );
    ST(2,:) = S0b*exp( (r-0.5*sigma_b^2)*T + sigma_b*sqrt(T).*Z(2,:) );
    ST(3,:) = S0c*exp( (r-0.5*sigma_c^2)*T + sigma_c*sqrt(T).*Z(3,:) );
    
    f_x = exp(-r*T)*max(sum(ST)- 200, 0);
    g_x = exp(-r*T)*max(ST(3,:)- S0c, 0);
    
    cov_fg = cov(f_x(1:n0),g_x(1:n0));
    alpha = -cov_fg(1,2)/cov_fg(2,2);
    
    b_option_control = -alpha*c0_c + mean(f_x(n0:end)+alpha.*g_x(n0:end));
    b_option = mean(f_x(n0:end));
    c_option = mean(g_x(n0:end));
    
    mat(idx,:) = [b_option_control, b_option, c_option ];
end

figure()
hold on
plot([0 50], [c0_c c0_c], 'k-') 
plot(mat(:,3), 'b*')
plot(mat(:,2),'bx')
plot(mat(:,1),'b.')
ylim([24 44])
legend({'Black Scholes Call', 'Crude Monte Carlo Call',...
    'Crude Monte Carlo Basket', 'Control Variates Basket'}, 'Location',...
    'east')
hold off 
