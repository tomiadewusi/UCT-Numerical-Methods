clear
clc
rng(0)
% Computing price of the option at inception
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing the intial values
S0 = [35; 65];
sigma1 = 0.25;
sigma2 = 0.2;
mu1 = 0.12;
mu2 = 0.15;
Sigma = [1 0.6; 0.6 1];
r = 0.1;
T = 1/12;
N = 300;
deltat = T/N;
K = 100;
L = chol(Sigma,'lower');
n = 50001;

% We now have a Normalised Sobol Sequence of length n
tic
ZS = norminv([Sobol(gendirnums(7,[1 3],16),n);...
    Sobol(gendirnums(11,[1 1 5],16),n)]);
toc
ZS(:,1) = [];
ZS = L*ZS;
ST = zeros(2,50000);
ST(1,:) = S0(1)*exp((r-0.5*sigma1^2)*T + sigma1*sqrt(T).*ZS(1,:));
ST(2,:) = S0(2)*exp((r-0.5*sigma2^2)*T + sigma2*sqrt(T).*ZS(2,:));

% Price of the option under the Risk Neutral measure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Price = mean(exp(-r*T)*max(sum(ST) - K, 0)); % 2.6806
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the intial delta values
indicator = sum(ST(:,1:5000)) > K;
h01 = mean(exp(-r*T).*(indicator).*ST(1,1:5000)./S0(1));  % 0.5718
h02 = mean(exp(-r*T).*(indicator).*ST(2,1:5000)./S0(2));  % 0.5691
b0 = Price - h01*S0(1) - h02*S0(2); % -54.3239
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I know that everything above this line is above board
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Now we calculate all the bi values to get the hedge portfolio profit and
% loss
% Intialise repeated values outside of the for loop

num = 100; % This is the number of PnL values
PnL = zeros(num,1);
bi = zeros(1,N); % hi for i = 0,1,2,...N -> There are N+1 columns
hi = zeros(2,N); % hi for i = 0,1,2,...N -> There are N+1 columns
bi(1) = b0;
ZS = ZS(:,1:5000);
ST = zeros(2,5000);

% Calculate the real world Stock price paths for both Stocks
Z = L*randn(2,N);
Si1 = S0(1)*...
    exp(cumsum((mu1-0.5*sigma1.^2)*deltat + sigma1*sqrt(deltat).*Z(1,:),2));
Si2 = S0(2)*...
    exp(cumsum((mu2-0.5*sigma2.^2)*deltat + sigma2*sqrt(deltat).*Z(2,:),2));

ZS1 = ZS(1,:);
ZS2 = ZS(2,:);
mat1 = exp(( r - 0.5*sigma1^2)*T+ sigma1*sqrt(T).*ZS1);
mat2 = exp(( r - 0.5*sigma2^2)*T+ sigma2*sqrt(T).*ZS2);
tic
for ix = 1:N
    ST1 = Si1(ix).*mat1;
    ST2 = Si2(ix).*mat2;
    
    indicator_discounted = exp(-r*T).*(ST1+ST2 > K);
    hi(:,ix) = [mean(indicator_discounted.*ST1./Si1(ix));...
        mean(indicator_discounted.*ST2./Si2(ix))];
    
end
hi = [[h01; h02] hi]; 
hi'
