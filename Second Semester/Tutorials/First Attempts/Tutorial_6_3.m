clear
clc
rng(0)

S0 = 100;
nu0 = 0.06;
theta = 0.06;
kappa = 9;
sigma = 0.5;
rho = -0.4;
T = 0.5;
r = 0.03;
K = 105;
Nm = 10;
Nq = 100;
umax = 30;


% Characteristic function pricing
charST =@(u) HestCharST(u,S0,nu0,theta,kappa,sigma,rho,T,r);

delta_u = umax/Nq;
un = ((1:Nq) - 0.5).*delta_u;
k = log(K);
P1 = 0.5 + (1/pi).*sum(    real(   (  exp(-1i.*un.*k).*charST(un-1i)  )...
    ./(  1i.*un.*charST(-1i)  )   ).*delta_u    );
P2 = 0.5 + (1/pi).*sum(real((exp(-1i.*un.*k).*charST(un)...
    ./(1i.*un))).*delta_u);
C_CF = S0*P1 - K*exp(-r*T)*P2;

P1_int = real(   (  exp(-1i.*un.*k).*charST(un-1i)  )...
    ./(  1i.*un.*charST(-1i)  )   )  ;
P2_int = real((exp(-1i.*un.*k).*charST(un)...
    ./(1i.*un))) ;

% Milstein Method Pricing
S = zeros(50,2);
delta_t = T/Nm;
for ix = 1:50
    n = 1000*ix;
    sOld = S0*ones(1,n);
    nuOld = nu0*ones(1,n);
    for jx = 1:Nm
        z = randn(2,n); % don't forget to generate this everytime 
        z1 = z(1,:);
        z2 = z(2,:);
        z2 = rho.*z1 + sqrt(1 - rho.^2).*z2;
        nuNew = nuOld + kappa.*(theta - nuOld).*delta_t...
            + sigma.*sqrt(nuOld).*sqrt(delta_t).*z2 ...
            + 0.25.*sigma.^2.*(z2.^2 -1).*delta_t;
        nuNew = max(nuNew,0);
        sNew = sOld.*exp((r-0.5.*nuOld).*delta_t...
            + sqrt(nuOld).*sqrt(delta_t).*z1);
        sOld = sNew;
        nuOld = nuNew;
    end
    sample = exp(-r*T).*max(sNew - K,0);
    S(ix,:) = [mean(sample) 3*std(sample)/sqrt(n)];
end

lower = C_CF - S(:,2);
upper = C_CF + S(:,2);

figure()
hold on
plot(lower, 'b--')
plot(upper,'b--')
plot(S(:,1),'b.')
plot([0 50], [C_CF C_CF],'k-')
hold off

figure()
hold on
plot(un,P1_int)
plot(un,P2_int)
hold off
