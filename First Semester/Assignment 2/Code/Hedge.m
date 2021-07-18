function [value1, value2] = Hedge
% Computing price of the option at inception
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing the intial values
S0 = [35; 65];
sigma = [0.25; 0.2];
mu = [0.12;0.15];
Sigma = [1 0.6; 0.6 1];
r = 0.1;
T = 1/12;
N = 300;
h = T/N;
K = 100;
L = chol(Sigma,'lower');
n = 50001;

% Sobol Sequence of length n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ZS = norminv([Sobol(gendirnums(7,[1 3],16),n);...
    Sobol(gendirnums(11,[1 1 5],16),n)]);

ZS = L*ZS(:,2:n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ST = zeros(2,50000);
ST(1,:) = S0(1)*exp((r-0.5*sigma(1)^2)*T + sigma(1)*sqrt(T).*ZS(1,:));
ST(2,:) = S0(2)*exp((r-0.5*sigma(2)^2)*T + sigma(2)*sqrt(T).*ZS(2,:));

% Price of the option under the Risk Neutral measure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Price = mean(exp(-r*T)*max(sum(ST) - K, 0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the intial delta values
indicator = exp(-r*T).*(sum(ST(:,1:5000)) > K);
h1_1 = mean(indicator.*ST(1,1:5000)./S0(1));
h1_2 = mean(indicator.*ST(2,1:5000)./S0(2));
b1 = Price - h1_1*S0(1) - h1_2*S0(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we calculate all the bi values to get the hedge portfolio profit and
% loss
% Intialise repeated values outside of the for loop

PnL = zeros(10000,1);
bi = zeros(1,N); % hi for i = 0,1,2,...N -> There are N+1 columns
hi = zeros(2,N); % hi for i = 0,1,2,...N -> There are N+1 columns
bi(1) = b1;
hi(:,1) = [h1_1; h1_2];
ZS = ZS(:,1:5000);
ST = zeros(2,5000);
var1 = (mu(1)-0.5*sigma(1).^2)*h;
var2 = sigma(1)*sqrt(h);
var3 = (mu(2)-0.5*sigma(2).^2)*h;
var4 = sigma(2)*sqrt(h);
var5 = (r-0.5*sigma(1)^2);
var6 = (r-0.5*sigma(2)^2);
var7 = exp(r*h);
ZS1 = ZS(1,:);
ZS2 = ZS(2,:);
for jx = 1:10000
    % Calculate the real world Stock price paths for both Stocks
    Z = L*randn(2,N);
    Si1 = S0(1)*exp(cumsum( var1 + var2.*Z(1,:),2));
    Si2 = S0(2)*exp(cumsum( var3 + var4.*Z(2,:),2));
    for ix = 2:N
        % Calculate the T value
        T = (1/12) - h*(ix - 1) ;
        ST(1,:) = Si1(ix-1)*exp(var5*T + sigma(1)*sqrt(T).*ZS1);
        ST(2,:) = Si2(ix-1)*exp(var6*T + sigma(2)*sqrt(T).*ZS2);
        indicator = exp(-r*T).*(sum(ST) > K);
        hi(:,ix) = [mean(indicator.*ST(1,:)./Si1(ix-1));...
            mean(indicator.*ST(2,:)./Si2(ix-1))];
        
        bi(ix) = bi(ix-1)*var7 - (hi(1,ix) - hi(1,ix-1))*Si1(ix-1) ...
            - (hi(2,ix) - hi(2,ix-1))*Si2(ix-1);
    end
    PnL(jx) = bi(N)*var7 + hi(1,N)*Si1(N) + hi(2,N)*Si2(N) -...
        max(Si1(N) + Si2(N) - K, 0);
end
value1 = Price;
value2 = PnL;
end