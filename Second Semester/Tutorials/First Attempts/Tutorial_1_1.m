clear 
clc 
rng(0) 

Nplus = 12; 
Nminus = -12; 
M = 20; 
T = 0.8; 
deltax = pi/Nplus; 
deltatau = T/M; 
theta = 0.5; 
% theta = 0; 
kappa = deltatau/(deltax.^2);
u0 = @(x) abs(sin(x)); 
um = @(x,t) t; 
up = @(x,t) t; 

dim = Nplus - Nminus -1; % This because we take the boundary conditions as 
                         % given 


C1 = (1 + 2*theta*kappa).*eye(dim,dim) ; 
C2 = diag((-theta*kappa).*ones(1,dim-1),1); 
C3 = diag((-theta*kappa).*ones(1,dim-1),-1);
C = C1 + C2 + C3; 

D1 = (1 - 2*(1-theta)*kappa).*eye(dim,dim) ; 
D2 = diag(((1-theta)*kappa).*ones(1,dim-1),1); 
D3 = diag(((1-theta)*kappa).*ones(1,dim-1),-1);
D = D1 + D2 + D3; 

U = zeros(dim,M+1) ; % We go from time 0 to M and MATLAB starts at 1
U(:,1) = u0(deltax.*((Nminus+1):(Nplus-1))).'; % Intial condition

B = [um(Nminus*deltax,deltatau*(0:M)); 
    zeros(dim-2,M+1); 
    up(Nplus*deltax,deltatau*(0:M))]; 

for ix = 2:M+1
    bm = B(:,ix-1); 
    bm_plus1 = B(:,ix); 
    U(:,ix) = C\(D*U(:,ix-1) + kappa*((1-theta)*bm + theta*bm_plus1));
    
end 

t = 0:deltatau:T; 
x = linspace(Nminus,Nplus,dim); 
[X,Tt] = meshgrid(t,x) ; 
figure()
surf(Tt,X,U)
