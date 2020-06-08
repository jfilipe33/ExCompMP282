function drdt = MAVtrans(t,m,D,g,,e2,r,f1,f2)
    drdt = zeros(2,1);
    drdt(1) = r(2);
    drdt(2) = -g*e2 + D*e2*(f1+f2)/m;
