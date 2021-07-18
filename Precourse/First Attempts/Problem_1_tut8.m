clear
clc

n = 10000;
rng(1, 'v4');
U = randn(2,n);

p = 0.75;

x1 = U(1,:);
x2 = p.*U(1,:) + sqrt(1-p.^2).*U(2,:);

X = [x1;x2];

%Principal component analysis
%Step 1 - mean

m = sum(X,2)./n;

%Step 2 - subtract away the mean 

m_subtract = repmat(m,1,n);
X = X - m_subtract;

%Step 3 - covariance 
%BEDMAS is always important
%always use brackets 

C = (1/(n-1)).*X*X.';

%Step 4 - eigenvectors and eigenvalues 

[W,D] = eig(C);

%Step 5 - orthogonal transform 

Y = W.'*X;

%First Plot  

figure()
subplot(1,2,1)

plot(X(1,:), X(2,:), 'ko')
axis equal
hold on 

%preparing the eigenvectors to be plotted 

%plot([0,D(1,1)*W(1,1)], [0,D(1,1)*W(2,1)])
w1 = [[0;0] D(1,1).*W(:, 1)]; 
w2 = [[0;0] D(2,2).*W(:, 2)]; 


plot(w1(1,:), w1(2,:), 'r-')
plot(w2(1,:), w2(2,:), 'r-')

axis equal
hold off
%Second plot
 
subplot(1,2,2)
plot(Y(1,:), Y(2,:), 'ko')
axis equal 
hold on 
%preparing the eigenvectors to be plotted 
% w1 = [[0;0] D(1,1).*W(:, 1)]; 
% w2 = [[0;0] D(2,2).*W(:, 2)]; 

% plot(w1(1,:), w1(1,:), 'r-')
% plot(w1(2,:), w1(2,:), 'r-')
hold off 



%note that axis equal must go before the hold off command 


%subplot(1,2,2)
% figure()
% subplot(1,2,1)
% plot(X(1,:), Y(1,:), 'bo')

