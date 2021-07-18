clear
clc
rng(0)

n = 50000;
index = 1:n;

sq3 = sqrt(3*(index - 1).^2) - floor(sqrt(3*(index - 1).^2));
sq7 = sqrt(7*(index - 1).^2) - floor(sqrt(7*(index - 1).^2));

nbins = 100;
figure()
subplot(1,2,1)
hist(sq3,nbins)
subplot(1,2,2)
hist(sq7,nbins)

n = 1000;
index = 1:n;
sq3 = sqrt(3*(index - 1).^2) - floor(sqrt(3*(index - 1).^2));
sq7 = sqrt(7*(index - 1).^2) - floor(sqrt(7*(index - 1).^2));
figure()
scatter(sq3, sq7,'b.')
xlabel('sqrt 3')
ylabel('sqrt 7')

estimates = zeros(50,1);
for idx = 1:50
    n = 1000*idx;
    index = 1:n;
    sq3 = sqrt(3*(index - 1).^2) - floor(sqrt(3*(index - 1).^2));
    sq7 = sqrt(7*(index - 1).^2) - floor(sqrt(7*(index - 1).^2));
    mat = zeros(10,10); 
    for ix = 1:10
        for jx = 1:10
            mat(ix,jx) = abs(mean((sq3 < jx*0.1).*(sq7 < ix.*0.1))...
                -(jx*0.1)*(ix*0.1));  
        end 
    end
    estimates(idx,1) = max(max(mat)); 
end
figure()
plot(1:50, estimates)