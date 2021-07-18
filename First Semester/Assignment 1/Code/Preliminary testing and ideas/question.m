clear
clc

mat1 = zeros(32,32);
mat2 = zeros(32,32);

% sigmaV = sort(0.1 + (0.3-0.1).*rand(1,32)); 
% KV = sort(80 + (120-80).*rand(1,32));
sigmaV = linspace(0.1,0.3,32); 
KV = linspace(80,120,32); 
eta = -1;
n = 15000;

for ix = 1:32
    for jx = 1:32
        mat1(ix,jx) = exact_soln(sigmaV(ix), KV(jx), eta);
    end
end

for ix = 1:32
    for jx = 1:32
        mat2(ix,jx) = control_anti_soln(sigmaV(ix), KV(jx), eta,n);
    end
end

error_mat1 = abs(mat1-mat2);
error_mat2 = abs(1 - mat2./mat1);
am_error = sum(sum(error_mat1,'omitnan'),'omitnan');
num_errors = sum(sum(error_mat1 > 0.05));

[X,Y] = meshgrid(KV,sigmaV);

figure()
surf(X,Y,error_mat1)
colormap(parula)  
xlabel('Strike Price')
ylabel('Volatility')
zlabel('Absolute value of error of Option')
title('Absolute Difference between my and exact solution: eta = -1')

% figure()
% surf(X,Y,error_mat2)
% xlabel('Strike Price')
% ylabel('Volatility')
% zlabel('Error as a percentage of Option')
% title('Percentage difference between my and exact solution: eta = -1')
% figure()
% surf(X,Y,mat1)
% xlabel('Strike Price')
% ylabel('Volatility')
% zlabel('Value of Option')
% title('Exact solution: eta = -1')

