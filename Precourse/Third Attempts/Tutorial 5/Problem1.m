clear 
clc 

%rng(1,'v4') not yet implemented in octave 
n = 10000; 
p = 0.75;  
U = randn(2,n); 
u1 = U(1,:); 
u2 = U(2,:); 

x1 = u1; 
x2 = p.*u1 + sqrt(1-p^2).*u2; 

X = [x1;x2]; 
m = mean(X,2); 

X = X - repmat(m,1,n); 

C = (1/(1-n)).*(X*X.'); 

[W,D] = eig(C);
w1 = [ [0;0] D(1,1).*W(:,1) ]; 
w2 = [ [0;0] D(2,2).*W(:,2) ]; 

Y = W.'*X;
C1 = (1/(1-n)).*(Y*Y.');
[W1,D1] = eig(C1);
w3 = [ [0;0] D1(1,1).*W1(:,1) ]; 
w4 = [ [0;0] D1(2,2).*W1(:,2) ]; 
 
figure()
subplot(1,2,1)
hold on 
plot(X(1,:), X(2,:), 'ko')
plot(w1(1,:),w1(2,:), 'r-')
plot(w2(1,:),w2(2,:), 'r-')
axis equal 
hold off 


subplot(1,2,2)
hold on 
plot(Y(1,:), Y(2,:), 'ko')
plot(w3(1,:),w3(2,:), 'b-')
plot(w4(1,:),w4(2,:), 'b-')
axis equal 
hold off 
