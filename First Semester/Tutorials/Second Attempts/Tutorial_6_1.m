clear
clc
rng(0)

xmax = 4;
xmin = 2;
ymax = 3;
ymin = 1;
mat1 = zeros(100,2); 
fx = @(U) (xmax-xmin).*(ymax-ymin).*sin((xmin + (xmax-xmin).*U(:,1))...
        .*(ymin + (ymax-ymin).*U(:,2)) );
for ix = 1:100
    n =1000*ix; 
    U = rand(n,2);
    fxvalues = fx(U); 
    mat1(ix,:) = [mean(fxvalues) 3*std(fxvalues)/sqrt(n)] ; 
end
fun = @(x,y) sin(x.*y);
q = integral2(fun,xmin,xmax,ymin,ymax); 

lower = q - mat1(:,2); 
upper = q + mat1(:,2); 

figure()
hold on 
plot([0 100], [q q], 'k-') 
plot(lower,'b--')
plot(upper,'b--') 
plot(mat1(:,1),'b.') 
hold off
