clear
clc
rng(0)

T = 1;
tdash = 0;
N = 6;
h = (T-tdash)/N;
sigma = 0.45;
r = 0.12;
S0 = 100;
K = 50;
ratio = 0.1;

mu_bar = (r - 0.5*sigma^2)*tdash + (r - 0.5*sigma^2)*(0.5*(T-tdash));
sigma_bar = sqrt(tdash*sigma^2+...
    ((T-tdash)/6)*(((T-tdash)/(T-tdash+h))+1)*sigma^2);

temp1 = S0*exp(mu_bar + 0.5*sigma_bar^2);
temp2 = normcdf( ((log(S0/K)+ mu_bar)/(sigma_bar)) +sigma_bar);
temp3 = K*normcdf( (log(S0/K) + mu_bar)/(sigma_bar));
c_geom = exp(-r*T)*(temp1*temp2 - temp3);

mat = zeros(50,3); 
if tdash == 0
    for idx = 1:50
        n = 1000*idx;
        n0 = ratio*n;
        Z = randn(N,n);
        Sj = [S0.*ones(1,n);S0*exp(cumsum((r-0.5*sigma.^2)...
            *h + sigma*sqrt(h).*Z))];
        f_x = exp(-r*T)*max(mean(Sj) - K, 0);
        g_x = exp(-r*T)*max(geomean(Sj) - K, 0);
        
        cov_fg = cov(f_x(1:n0), g_x(1:n0));
        alpha = -cov_fg(1,2)/cov_fg(2,2);
        
        monte_geo = mean(f_x(n0+1:end)); 
        monte_arth = mean(g_x(n0+1:end));
        control_c_arth = monte_geo...
            + alpha*(monte_arth - c_geom);
        
        mat(idx,:) = [monte_geo monte_arth control_c_arth];
    end
else
    for idx = 1:50
        n = idx*1000;
        n0 = ratio*n;
        Z =  randn(N+1,n);
        
        j = repmat([0:N]',1,n);
        Z0 = repmat(Z(:,1),1,n);
        Zk = cumsum([zeros(N+1,1) Z(:,2:end)]);
        
        Sj = S0*exp((r-0.5*sigma^2)*(tdash + j.*h)+...
            sigma*(sqrt(tdash).*Z0+sqrt(h).*Zk));
        f_x = exp(-r*T)*max(mean(Sj) - K, 0);
        g_x = exp(-r*T)*max(geomean(Sj) - K, 0);
        
        cov_fg = cov(f_x(1:n0), g_x(1:n0));
        alpha = -cov_fg(1,2)/cov_fg(2,2);
        
        monte_geo = mean(f_x(n0+1:end)); 
        monte_arth = mean(g_x(n0+1:end));
        control_c_arth = monte_geo...
            + alpha*(monte_arth - c_geom);
        mat(idx,:) = [monte_arth monte_geo control_c_arth];
    end
end

figure()
hold on 
plot([0 50], [c_geom c_geom])
plot(mat(:,1),'bx')
plot(mat(:,2),'b*')
plot(mat(:,3),'b.')
ylim([47 50.5])
hold off
