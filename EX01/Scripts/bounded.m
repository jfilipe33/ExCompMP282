function b = bounded(rx,ry,wx,wy)
    taux=0:.01:2*pi;
    x=cos(taux)*0.1+wx;
    y=sin(taux)*0.1+wy;
    b = inpolygon(rx,ry,x,y);