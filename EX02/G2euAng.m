function ang = G2euAng(obj)
      
    obj.DTilde =  ((1 - obj.gTilde'*obj.gTilde)*eye(3) + 2*obj.gTilde*obj.gTilde' - 2*skew(obj.gTilde))/(1 + obj.gTilde'*obj.gTilde);
    obj.D = (inv(obj.Dbar)*obj.DTilde)';
    phi = -atand(obj.D(3,2)/obj.D(3,3));
    theta = asind(obj.D(3,1));
    psi = -atand(obj.D(2,1)/obj.D(1,1));
    ang = [phi theta psi];
    %obj.Omegabar = [cosd(theta)*cosd(psi) sind(psi) 0; -cosd(theta)*sind(psi) cosd(psi) 0; sind(theta) 0 1];
        
end
