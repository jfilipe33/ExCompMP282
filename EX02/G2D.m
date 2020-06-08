function obj = G2D(obj)
      
    obj.DTilde =  ((1 - obj.gTilde'*obj.gTilde)*eye(3) + 2*obj.gTilde*obj.gTilde' - 2*skew(obj.gTilde))/(1 + obj.gTilde'*obj.gTilde);
    obj.D = (inv(obj.Dbar)*obj.DTilde)';
        
end
