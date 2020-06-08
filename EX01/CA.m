function f1f2 = CA(f_c,tau_c,l)
    f1f2 = zeros(2,1);
    f1f2 = [0.5 1/(2*l);0.5 -1/(2*l)]*[norm(f_c);tau_c];