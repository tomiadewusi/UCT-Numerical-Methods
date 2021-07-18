function W = BrownianBridge(t,W0,WN,Z)
	if isempty(Z) 
		W = [];
    else
	
        N = length(t)-1;
        m = floor((N+1)*0.5);

        Wm = (((t(end)-t(m+1))*W0+(t(m+1)-t(1))*WN)/(t(end)-t(1)))...
    + sqrt((t(m+1)-t(1))*(t(end)-t(m+1))/(t(end)-t(1)))*Z(m);

        Wl = BrownianBridge(t(1:m+1),W0,Wm,Z(1:m-1));
        Wr = BrownianBridge(t(m+1:end),Wm,WN,Z(m+1:end));
        W = [Wl Wm Wr]; 
    end 
end