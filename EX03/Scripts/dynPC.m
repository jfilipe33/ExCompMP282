function states = dynPC(st,obj) 
    
    r = st(1:3);
    v = st(4:6);
    e3 = [0 0 1]';
    g = 9.81;
    m = 0.5;

    states = [v;
              (1/m)*obj.Fc - g*e3];
          
end