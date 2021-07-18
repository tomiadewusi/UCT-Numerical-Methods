clear 
clc
rng(0) 

sigma = 0.4; 
r = 0.06; 
K = 50; 
T = 1; 

Nminus = -100; 
Nplus = 20; 
M = 35; 
deltax = 0.06; 
deltatau = (T*sigma.^2)/(2*M) ; 
theta = 0.5; 
kappa = deltatau/(deltax.^2);

gamma = 2*r/(sigma.^2); 
beta = -0.25*(gamma + 1).^2; 
alpha = -0.5*(gamma - 1); 


fs = @(x) K.*exp(x); 
ft = @(t) T - 2.*t./(sigma.^2); 
f = @(V,x,t) (V./K).*exp(-alpha.*x - beta.*t); 

VT = @(S) (S > K); 
V0 = @(S,t) 0; 
Vinf = @(S,t) exp(-r.*(T-t)); 

u0 = @(x) f(VT(fs(x)),x,0); 
um = @(x,t) f( V0(fs(x),ft(t)) ,x ,t ); 
up = @(x,t) f( Vinf(fs(x),ft(t)) ,x ,t );

dim = Nplus - Nminus -1; 

C1 = (1 + 2*theta*kappa).*eye(dim,dim) ; 
C2 = diag((-theta*kappa).*ones(1,dim-1),1); 
C3 = diag((-theta*kappa).*ones(1,dim-1),-1);
C = C1 + C2 + C3; 

D1 = (1 - 2*(1-theta)*kappa).*eye(dim,dim) ; 
D2 = diag(((1-theta)*kappa).*ones(1,dim-1),1); 
D3 = diag(((1-theta)*kappa).*ones(1,dim-1),-1);
D = D1 + D2 + D3; 

U = zeros(dim,M+1) ; 
U(:,1) = u0(deltax.*linspace(Nminus+1,Nplus-1,dim)).'; 

B = [um(Nminus*deltax,deltatau*(0:M)); 
    zeros(dim-2,M+1); 
    up(Nplus*deltax,deltatau*(0:M))]; 

for ix = 2:M+1
    bm = B(:,ix-1); 
    bm_plus1 = B(:,ix); 
    U(:,ix) = C\(D*U(:,ix-1) + kappa*((1-theta)*bm + theta*bm_plus1));   
end 
% You include intial conditions but you exclude boundary conditions 
t = linspace(0,M,M+1); 
x = linspace(Nminus+1,Nplus-1,dim); 
[Tt,X] = meshgrid(t,x) ; 
figure()
surf(X,Tt,U)
zlabel('X')
title('diffusion solution') 

Sn = fs(x.*deltax); 
tm = ft(t.*deltatau); 
[M,N ]= meshgrid(tm,Sn) ; 

V = K.*exp(alpha.*X*deltax + beta.*Tt.*deltatau).*U ; 
figure()
surf(M,N,V)
xlabel('tau')
ylabel('S')
zlabel('V')

title('call solution')


S = Sn; 
d = (log(S/K) + (r - 0.5.*sigma.^2).*T)./(sigma.*sqrt(T)); 
exact = exp(-r*T)*normcdf(d); 

diffusion = (V(:,end)).'; 
error = exact - diffusion ; 
figure()
plot(S,error) 