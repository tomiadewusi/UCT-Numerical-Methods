function vector = virus_strain1(sigma, mus, mu1, mu2, k2, t, v)
    temp = zeros(3,1);
    
    if t < 365
        k1 = 1/2; 
    else 
        k1 = 1/3; 
    end 
    
    temp(1) = sigma - mus.*v(1) - v(1)*( k1*v(2) + k2*v(3) ); 
    temp(2) = k1*v(1)*v(2) - mu1*v(2); 
    temp(3) = k2*v(1)*v(3) - mu2*v(3); 
    
    vector = temp; 