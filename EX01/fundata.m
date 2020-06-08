function [Os,tp,ts] = fundata(k1,k2)
    wn = sqrt(k1);
    eta = k2/(2*wn);
    Os = exp((-pi*eta)/sqrt(1-eta^2))
    tp = pi/(sqrt(1-eta^2)*wn)
    ts = -log(0.05)/(eta*wn)