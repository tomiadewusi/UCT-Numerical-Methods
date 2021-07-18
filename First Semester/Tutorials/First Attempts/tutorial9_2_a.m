clear
clc
rng(0)

n = 3000;
Z = randn(n,1);
nbins = 50;
hist(Z,nbins)

% r = a + (b-a).*rand(N,1).
d = 50;

X = (((1:50) - 1)/d) + (1/d).*rand(n/d,d);

X = norminv(reshape(X,[],1));
figure()
subplot(1,2,1)
hist(Z,nbins)
subplot(1,2,2)
hist(X,nbins)
axis tight
