function obj = DRefresh(obj)
    phi = obj.alpha(1);
    theta = obj.alpha(2);
    psi  = obj.alpha(3);
    
    obj.D = [cosd(psi)*cosd(theta) (cosd(psi)*sind(theta)*sind(phi))+(sind(psi)*cosd(phi)) -(cosd(psi)*sind(theta)*cosd(phi))+(sind(psi)*sind(phi));
        -sind(psi)*cosd(theta) -(sind(psi)*sind(theta)*sind(phi))+(cosd(psi)*cosd(phi)) (sind(psi)*sind(theta)*cosd(phi))+(cosd(psi)*sind(phi));
        sind(theta) -cosd(theta)*sind(phi) cosd(theta)*cosd(phi)];
    
        
end