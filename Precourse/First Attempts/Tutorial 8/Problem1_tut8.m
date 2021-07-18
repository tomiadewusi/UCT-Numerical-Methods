clear 
clc

f = @(x) (1/(sqrt(2.*pi)) ) .* exp( -((x.^2)/(2)) ); 
a = -3; 
b = 3; 
h = 0.1; 

format long 
mpoint = midpoint_tut8(f,a,b-h,h)
comp = composite_tut8(f,a,b,h)
checkQ = quad(f,a,b)
checkI = integral(f,a,b)


points = a:h:b; 
f_grad = gradient(f(points), h);
f_forward = f_dash_forward(f,a,b,h);
f_central = f_dash_central(f,a,b,h);

figure()
hold on 

plot(points,f(points))
plot(points,f_grad)
plot(points, f_forward)
plot(points, f_central)
legend('f',"f_g",'f_f','f_c')

hold off 

