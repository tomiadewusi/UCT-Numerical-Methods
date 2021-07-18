clear
clc

format long
N = 25;
Z0 = rand(1,N+1);
% Remember that you can handle function defaults by using anonmous
% functions
handle = @fun;
f = @(Z) fun(Z,K,sigma,eta,S0,L,U,r,T,deltat); 
v = fminunc(handle,Z0);
v = v/norm(v);
v = v.';
K1 = 100000;
final = zeros(K1,N+1); 
for idx = 1:K1
    U1 = rand(1);
    V = (idx - 1 + U1)/K1;
    X = norminv(V);
    Z = randn(N+1,1);
    I = eye(N+1,N+1);
    final(idx,:) = I*v.*X + (I - I*(v*v.')*I)*Z; 
end
K = 80; sigma = 0.2677; eta = 1; n = 100000; 
S0 = 100; L = 70; U = 130; r = 0.1; T = 1; N = 25; deltat = T/(N+1);
Z = final; 
Si = S0*exp(cumsum((r-0.5*sigma^2)*deltat + sigma*sqrt(deltat).*Z,2));
ST = Si(:,end );
Si = Si(:,1:end-1);
approx = mean(exp(-r*T).*(max(eta.*(ST - K),0)).*mean(Si > L & Si < U,2));

n1 = 100000; 
crude = crude_soln(sigma,K,eta,n1);

exact = exact_soln(sigma,K,eta);

control = control_anti_soln(sigma,K,eta,n); 

im_diff = abs(approx - exact)
control_diff = abs(control - exact) 
crude_diff = abs(crude - exact)
