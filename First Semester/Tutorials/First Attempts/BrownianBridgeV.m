function W = BrownianBridgeV(t,W0,WN,Z)

if isempty(Z) 
    W = [];
else

    N = length(t)-1;
    m = floor((N+1)*0.5);

    Wm = (((t(1,end)-t(1,m+1))*W0+(t(1,m+1)-t(1,1))*WN)/(t(1,end)-t(1,1)))...
       + sqrt((t(1,m+1)-t(1,1))*(t(1,end)-t(1,m+1))/(t(1,end)-t(1,1))).*Z(:,m);

    Wl = BrownianBridgeV(t(1:m+1),W0,Wm,Z(:,1:m-1));
    Wr = BrownianBridgeV(t(m+1:end),Wm,WN,Z(:,m+1:end));
    W = [Wl Wm Wr]; 
end 
end