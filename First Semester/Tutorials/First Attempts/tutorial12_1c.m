clear
clc
rng(0)

sigma_mat = [1, 0.8,0.9;
            0.8, 1, 0.75;
            0.9 0.75, 1];
L = chol(sigma_mat, 'lower');

S0a = 25;
sigma_a = 0.4;

S0b = 65;
sigma_b = 0.25;

S0c = 120;
sigma_c = 0.3;

K = 210;

T = 1.5;
r = 0.1;

% Generating the closed sequence
Z2 = vanderCorput(50001,2)';
Z3 = vanderCorput(50001,3)';
Z5 = vanderCorput(50001,5)';
hal = [Z2(2:end); Z3(2:end); Z5(2:end)];

mat = zeros(50,2);
mat2 = zeros(50,2);
for idx = 1:50
    n = 1000*idx;
    % Crude Monte Carlo
    Z = L*randn(3,n);
    
    ST = zeros(3,n);
    ST(1,:) = S0a*exp( (r-0.5*sigma_a^2)*T + sigma_a*sqrt(T).*Z(1,:) );
    ST(2,:) = S0b*exp( (r-0.5*sigma_b^2)*T + sigma_b*sqrt(T).*Z(2,:) );
    ST(3,:) = S0c*exp( (r-0.5*sigma_c^2)*T + sigma_c*sqrt(T).*Z(3,:) );
    
    f_x = exp(-r*T)*max(K-sum(ST), 0);
    option_value = mean(f_x);
    confidence = 3*std(f_x)/sqrt(n);
    mat(idx,:) = [option_value confidence];
    
    % Open Rule
    Z = L*norminv(hal(:,1:n));
    ST = zeros(3,n);
    ST(1,:) = S0a*exp( (r-0.5*sigma_a^2)*T + sigma_a*sqrt(T).*Z(1,:) );
    ST(2,:) = S0b*exp( (r-0.5*sigma_b^2)*T + sigma_b*sqrt(T).*Z(2,:) );
    ST(3,:) = S0c*exp( (r-0.5*sigma_c^2)*T + sigma_c*sqrt(T).*Z(3,:) );
    
    f_x = exp(-r*T)*max(K-sum(ST), 0);
    option_value = mean(f_x);
    mat2(idx,1) = option_value;
    
    % Closed Rule
    index = (1:n)./(n+1);
    ham = [index; Z2(2:n+1); Z3(2:n+1)];
    Z = L*norminv(ham);
    ST = zeros(3,n);
    ST(1,:) = S0a*exp( (r-0.5*sigma_a^2)*T + sigma_a*sqrt(T).*Z(1,:) );
    ST(2,:) = S0b*exp( (r-0.5*sigma_b^2)*T + sigma_b*sqrt(T).*Z(2,:) );
    ST(3,:) = S0c*exp( (r-0.5*sigma_c^2)*T + sigma_c*sqrt(T).*Z(3,:) );
    
    f_x = exp(-r*T)*max(K-sum(ST), 0);
    option_value = mean(f_x);
    mat2(idx,2) = option_value;
    
    
end

c = mat2(end,2);
lower = c - mat(:,2);
upper = c + mat(:,2);

figure()
hold on
plot([ 0 50], [c c], 'b-')
plot(mat(:,1),'b.')
plot([lower upper], 'b--')
plot(mat2(:,1),'r.')
plot(mat2(:,2),'k.')
hold off





