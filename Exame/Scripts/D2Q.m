function quat = D2Q(D)

    eta = 0.5*sqrt(1+trace(D));
    eps = 1/(4*eta)*[D(2,3)-D(3,2);D(3,1)-D(1,3);D(1,2)-D(2,1)];
    quat = [eps' eta]'; 
    
end