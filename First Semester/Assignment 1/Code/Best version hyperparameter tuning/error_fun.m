function [num_errors,am_error, tot_time,benchmark] ...
    = error_fun(testfun,eta,plotting,dim)
mat1 = zeros(dim,dim); 
mat2 = zeros(dim,dim);
mat3 = zeros(dim,dim); 
mat4 = zeros(dim,dim);  
sigmaV = linspace(0.1, 0.3, dim) ;
KV = linspace(80, 120, dim);

for ix = 1:dim
    for jx = 1:dim
        t1 = tic; 
        mat1(ix,jx) = exact(sigmaV(ix), KV(jx), eta);
        mat3(ix,jx) = toc(t1); 
    end
end

for ix = 1:dim
    for jx = 1:dim
        t2 = tic; 
        mat2(ix,jx) = testfun(sigmaV(ix), KV(jx), eta);
        mat4(ix,jx) = toc(t2); 
    end
end

error_mat = abs(mat1-mat2); 
am_error = sum(sum(error_mat,'omitnan'),'omitnan'); 
num_errors = sum(sum(error_mat > 0.05)); 
tot_time = sum(sum(mat4)); 
benchmark = tot_time/sum(sum(mat3)); 

if plotting == 1
    [X,Y] = meshgrid(KV,sigmaV);
    figure()
    surf(X,Y,error_mat)
    xlabel('Strike Price')
    ylabel('Volatility')
    zlabel('Absolute value of error')
end 
end

function value = exact(sigma, K, eta)
% Fixed Input values
S0 = 100; L = 70; U = 130; r = 0.1; T = 1; N = 25; deltat = T/(N+1);
t = deltat.*(1:25);

d1 = (log(S0/K) + (r+0.5*sigma^2)*T )/(sigma*sqrt(T)); 
d2 = d1 - sigma*sqrt(T); 

d3 = (log(S0/L) + (r+0.5*sigma^2).*t )./(sigma.*sqrt(t)); 
d4 = d3 - sigma.*sqrt(t); 

d5 = (log(S0/U) + (r+0.5*sigma^2).*t )./(sigma.*sqrt(t)); 
d6 = d5 - sigma.*sqrt(t); 

M = @(x,y,rho) mvncdf([x y],0,[1 rho; rho 1]);
p = -eta.*sqrt(t)./sqrt(T); 

tot = 0; 
for ix = 1:N 
    tot = tot + eta*(S0*(M(-d5(ix),eta*d1,p(ix))-M(-d3(ix), eta*d1,p(ix)))...
    -K*exp(-r*T)*(M(-d6(ix),eta*d2,p(ix))-M(-d4(ix), eta*d2,p(ix)))); 
end 
value = tot/N; 
end