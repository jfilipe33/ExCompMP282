function states = dynFull(st,obja,objp) 
    
    r_til = st(1:3);
    v_til = st(4:6);
    g_Til = st(7:9);
    W_Til = st(10:12);
    e3 = [0 0 1]';
    g = 9.81;
    m = 0.5;

    
    states = [v_til;
              -(1/m)*objp.Fc + g*e3 - (1/m)*objp.Fd + objp.aTar;
              0.5*(g_Til*g_Til' + skew(g_Til) + eye(3))*W_Til;
              -inv(obja.Jb)*skew(obja.Jb*W_Til)*W_Til - inv(obja.Jb)*(obja.Tc + obja.Td)];
          

end