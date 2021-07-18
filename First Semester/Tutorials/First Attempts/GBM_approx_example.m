clear 
clc 
rng(0)

mu = 0.2; 
sigma = 0.15; 
S0 = 100; 
t = 0:0.1:5; 
deltat = 0.1; 
Z = randn(length(t)-1, 1); 
%path
S_path =[S0;S0*exp(cumsum((mu-0.5*sigma^2)*deltat+sigma*sqrt(deltat).*Z ))]; 

%approximations 
S_euler = zeros(length(t),1); 
S_euler(1,:) = S0;
S_mils = zeros(length(t),1); 
S_mils(1,:) = S0; 

for idx = 2:length(t)
	S_euler(idx,:) = S_euler(idx-1,:)+mu.*S_euler(idx-1,:).*deltat...
		+sigma.*S_euler(idx-1,:).*sqrt(deltat).*Z(idx-1,:); 
end 

for idx = 2:length(t)
	S_mils(idx,:) = S_mils(idx-1,:)+mu.*S_mils(idx-1,:).*deltat...
		+sigma.*S_mils(idx-1,:).*sqrt(deltat).*Z(idx-1,:)...
		+0.5*sigma^2.*S_mils(idx-1,:).*(Z(idx-1,:)^2 - 1).*deltat; 
end 

%Plotting

figure()
subplot(1,2,1)
hold on 
plot(t,S_path)
plot(t,S_euler','linestyle','-','marker','.','markersize',12)
legend({'Exact solution','Euler-Maruyama approximation'},...
    'Location','northwest','NumColumns',1)
ylabel('S_t')
xlabel('t (years)')
axis tight
hold off 
subplot(1,2,2)
hold on 
plot(t,S_path)
plot(t,S_mils,'linestyle','-','marker','.','markersize',12)
legend({'Exact solution','Euler-Maruyama approximation'},...
    'Location','northwest','NumColumns',1)
ylabel('S_t')
xlabel('t (years)')
axis tight
hold off 
