clear 
clc 

randnums_1 = LCG(16807,0,(2^31)-1,1,10000); 
randnums_2 = LCG((2^16)+3,0,(2^31),1,10000); 

figure() 
subplot(1,2,1) 
hist(randnums_1) 
subplot(1,2,2) 
hist(randnums_2) 

u1_1 = randnums_1(1:2:end); 
u1_2 = randnums_1(2:2:end); 

u2_1 = randnums_2(1:2:end); 
u2_2 = randnums_2(2:2:end); 

figure()
subplot(1,2,1) 
scatter(u1_1,u1_2,'b.')
subplot(1,2,2) 
scatter(u1_1,u1_2,'b.')

u1_3d_1 = randnums_1(1:3:end-1); 
u1_3d_2 = randnums_1(2:3:end-1);
u1_3d_3 = randnums_1(3:3:end-1); 

u2_3d_1 = randnums_2(1:3:end-1); 
u2_3d_2 = randnums_2(2:3:end-1);
u2_3d_3 = randnums_2(3:3:end-1); 

figure() 
subplot(1,2,1) 
scatter3(u1_3d_1, u1_3d_2, u1_3d_3,'b.')
subplot(1,2,2) 
scatter3(u2_3d_1, u2_3d_2, u2_3d_3,'b.') 

