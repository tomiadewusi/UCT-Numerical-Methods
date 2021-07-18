clear 
clc
rng(0)

Z2 = vanderCorput(50001,2);
Z2 = Z2(2:end); 
Z3 = vanderCorput(50001,3);
Z3 = Z3(2:end); 
Z5 = vanderCorput(50001,5);
Z5 = Z5(2:end); 
hal = [Z2; Z3; Z5]; 

figure()
hold on 


subplot(2,3,1)
scatter(Z2(1:1000),Z3(1:1000),'b.')
subplot(2,3,4)
scatter(norminv(Z2(1:1000)),norminv(Z3(1:1000)),'b.')
title('Dimension 1 and 2')
 

subplot(2,3,2)
scatter(Z2(1:1000),Z5(1:1000),'b.')
subplot(2,3,5)
scatter(norminv(Z2(1:1000)),norminv(Z5(1:1000)),'b.')
title('Dimension 1 and 3')

subplot(2,3,3)
scatter(Z3(1:1000),Z5(1:1000),'b.')
subplot(2,3,6)
scatter(norminv(Z3(1:1000)),norminv(Z5(1:1000)),'b.')
title('Dimension 2 and 3')


hold off 


n = 1000; 
Z2 = vanderCorput(n+1,2)'; 
Z2 = Z2(2:end);
Z3 = vanderCorput(n+1,3)'; 
Z3 = Z3(2:end);

index = (1:n)./(n+1);
ham = [index; Z2; Z3];


figure()
subplot(2,3,1)
scatter(ham(1,:),ham(2,:),'b.')
subplot(2,3,4)
scatter(norminv(ham(1,:)),norminv(ham(2,:)),'b.')
title('Dimension 1 and 2')

subplot(2,3,2)
scatter(ham(1,:),ham(3,:),'b.')
subplot(2,3,5)
scatter(norminv(ham(1,:)),norminv(ham(3,:)),'b.')
title('Dimension 1 and 3')


subplot(2,3,3)
scatter(ham(2,:),ham(3,:),'b.')
subplot(2,3,6)
scatter(norminv(ham(2,:)),norminv(ham(3,:)),'b.')
title('Dimension 2 and 3')


hold off 

