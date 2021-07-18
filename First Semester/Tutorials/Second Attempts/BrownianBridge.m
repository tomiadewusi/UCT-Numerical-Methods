function W = BrownianBridge(t,W0,WN,Z)

if isempty(Z)
    W = []; 
else
    N = length(t) - 1; 
    m = floor((N+1)/2); 
    sum1 = ((t(N+1) - t(m+1)).*W0 + (t(m+1)-t(1)).*WN) / (t(N+1) - t(1)); 
    temp = (t(m+1)-t(1))*(t(N+1)-t(m+1))/(t(N+1)-t(1)); 
    sum2 = sqrt(temp).*Z(:,m);
    Wm = sum1 + sum2; 
    
    Wl = BrownianBridge(t(1:m+1), W0, Wm, Z(:,1:m-1)); 
    Wr = BrownianBridge(t(m+1:N+1), Wm, WN, Z(:,m+1:N-1)); 
    W = [Wl Wm Wr];
end 