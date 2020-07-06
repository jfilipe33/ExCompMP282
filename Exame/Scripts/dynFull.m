function states = dynFull(st,obja,objp,objm) 
    
    r = st(1:3);
    v = st(4:6);
    q = st(7:10);
    W = st(11:13);
    w = st(14:19);

    D = Q2D(st(7:10));

    
    states = [v;
              (1/objm.m)*D'*objm.e3*objp.Fc - objm.g*objm.e3 + (1/objm.m)*objp.Fd;
              0.5*matW(W)*q;
              inv(obja.Jb)*skew(obja.Jb*W)*W + inv(obja.Jb)*(obja.Tc + obja.Td);
              -1/objm.Tm*w+objm.km/objm.Tm*objm.wbar];         

end