function states = dynFull(st,obja,objp) 
    
    r = st(1:3);
    v = st(4:6);
    a = st(7:9);
    W = st(10:12);
    e3 = [0 0 1]';
    g = 9.81;
    m = 0.5;

    
    states = [v;
              (1/m)*objp.Fc - g*e3;
              matA(a)*W;
              inv(obja.Jb)*skew(obja.Jb*W)*W + inv(obja.Jb)*(obja.Tc)];
          

end