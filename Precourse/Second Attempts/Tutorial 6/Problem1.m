clear 
clc 
rng(1,'v4')

%Preparation 
n = 10000; 
U = randn(2,n); 
p = 0.75; 
u1 = U(1,:); 
u2 = U(2,:); 
x1 = u1; 
x2 = p.*u1 + sqrt(1 - p.^2).*u2; 

X = [x1;x2]; 

%Step 1 
m = mean(X,2); 

%Step 2 - removing the residual mean
X = X - repmat(m,1,n); 

%Step 3 - covariance 
C = (1./((n-1))).*X*X'; 

%Step 4 - eigenvectors and eigenvalues 
[W,D] = eig(C); 

e1 = [[0;0] D(1,1).*W(:,1)]; 
e2 = [[0;0] D(2,2).*W(:,2)]; 

%Step 5 - Orthogonal Transform 
Y = W'*X; 

Y(:,1)


%Plotting 
figure()
subplot(1,2,1)
plot(x1,x2,'ko')
hold on 
plot(e1(1,:),e1(2,:),'r-')
plot(e2(1,:),e2(2,:),'r-')
hold off 
axis equal % to ensure that you have a fixed aspect ratio 

subplot(1,2,2)
plot(Y(1,:),Y(2,:),'ko')
axis equal 
