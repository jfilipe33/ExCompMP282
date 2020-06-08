function xnew = dyn(x,mav)
                   
    v = x(3:4);
    W = x(6);
    g = 9.81;
    e2 = [0;1];

    xnew =  [v;
            (1/mav.m)*mav.D*e2*(mav.f(1)+mav.f(2)) - e2*g;
            W;
            (1/mav.J) * mav.Tc];
         
    mav.a = xnew(3:4);
    mav.Wdot = xnew(6);

end