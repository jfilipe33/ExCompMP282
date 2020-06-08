function states = dynAC(st,obj) 
    
    alpha = st(1:3);
    W = st(4:6);
    
    states = [matA(alpha)*W;
              inv(obj.Jb)*skew(obj.Jb*W)*W + inv(obj.Jb)*(obj.Tc)];

end