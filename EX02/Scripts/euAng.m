function obj = euAng(obj)
    phi = obj.angTar(1);
    theta = obj.angTar(2);
    psi  = obj.angTar(3);
    
    obj.Dbar = [cosd(psi)*cosd(theta) (cosd(psi)*sind(theta)*sind(phi))+(sind(psi)*cosd(phi)) -(cosd(psi)*sind(theta)*cosd(phi))+(sind(psi)*sind(phi));
        -sind(psi)*cosd(theta) -(sind(psi)*sind(theta)*sind(phi))+(cosd(psi)*cosd(phi)) (sind(psi)*sind(theta)*cosd(phi))+(cosd(psi)*sind(phi));
        sind(theta) -cosd(theta)*sind(phi) cosd(theta)*cosd(phi)];
    
    %obj.Omegabar = [cosd(theta)*cosd(psi) sind(psi) 0; -cosd(theta)*sind(psi) cosd(psi) 0; sind(theta) 0 1];
        
end
