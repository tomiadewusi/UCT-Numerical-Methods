function [value1,value2] = Hedge
% Computing price of the option at inception
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing the intial values
one = 1;
two = 2; 
zero = 0; 
fivethou = 5000; 
S01 =35; S02 = 65;
sigma1 = 0.25;
sigma2 = 0.2;
mu1 = 0.12;
mu2 = 0.15;
Sigma = [one 0.6; 0.6 one];
r = 0.1;
T = one/12;
T_index = T; 
N = 300;
h = T/N;
K = 100;
L = chol(Sigma,'lower');
n = 50001;
k = 16; 

% Sobol Sequence of length n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ZS = norminv([Sobol(gendirnums(7,[one 3],k),n);...
              Sobol(gendirnums(11,[one one 5],k),n)]);
ZS = L*ZS(:,two:n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ST1 = S01*exp((r-0.5*sigma1^2)*T + sigma1*sqrt(T).*ZS(one,:));
ST2 = S02*exp((r-0.5*sigma2^2)*T + sigma2*sqrt(T).*ZS(two,:));
% Price of the option under the Risk Neutral measure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Price = sum(exp(-r*T)*max(ST1 + ST2 - K, zero))/(n-one);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the intial delta values
ST1 = ST1(one:fivethou); 
ST2 = ST2(one:fivethou); 
indicator = exp(-r*T).*(ST1+ ST2 > K);
h1_1 = sum(indicator.*ST1./S01)/fivethou;
h1_2 = sum(indicator.*ST2./S02)/fivethou;
b1 = Price - h1_1*S01 - h1_2*S02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we calculate all the bi values to get the hedge portfolio profit and
% loss
%
% Intialise repeated values outside of the for loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num = 10; 
PnL = zeros(num,1);
bi = zeros(one,N); % hi for i = 0,1,2,...N -> There are N+1 columns
hi = zeros(two,N); % hi for i = 0,1,2,...N -> There are N+1 columns
bi(1) = b1;
hi(:,1) = [h1_1; h1_2];
var1 = (mu1-0.5*sigma1.^2)*h;
var2 = sigma1*sqrt(h);
var3 = (mu2-0.5*sigma2.^2)*h;
var4 = sigma2*sqrt(h);
var5 = (r-0.5*sigma1^2);
var6 = (r-0.5*sigma2^2);
var7 = exp(r*h);
ZS = ZS(:,one:fivethou); 
ZS1 = ZS(one,:).*sigma1; 
ZS2 = ZS(two,:).*sigma2;  
for jx = one:num
    % Calculate the real world Stock price paths for both Stocks
    Z = L*randn(two,N);
    Si1 = S01*exp(cumsum(var1 + var2.*Z(one,:),two));
    Si2 = S02*exp(cumsum(var3 + var4.*Z(two,:),two));
    for ix = two:N
        T = T_index - h*(ix - one) ;
        ST1 = Si1(ix-one)*exp(var5*T + sqrt(T).*ZS1);
        ST2 = Si2(ix-one)*exp(var6*T + sqrt(T).*ZS2);
        indicator = exp(-r*T).*(ST1 + ST2 > K);
        
        hi(:,ix) = [sum(indicator.*ST1./Si1(ix-one))/fivethou;...
                    sum(indicator.*ST2./Si2(ix-one))/fivethou];
        
        bi(ix) = bi(ix-one)*var7...
            - (hi(one,ix) - hi(one,ix-one))*Si1(ix-one) ...
            - (hi(two,ix) - hi(two,ix-one))*Si2(ix-one);
    end
    PnL(jx) = bi(N)*var7 + hi(1,N)*Si1(N) + hi(two,N)*Si2(N) -...
        max(Si1(N) + Si2(N) - K, zero);
end
value1 = Price;
value2 = PnL;
end