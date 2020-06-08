function obj = aTildeRefresh(obj)
    phi = obj.alpha(1);
    theta = obj.alpha(2);
    psi  = obj.alpha(3);
    
    D = [cosd(psi)*cosd(theta) (cosd(psi)*sind(theta)*sind(phi))+(sind(psi)*cosd(phi)) -(cosd(psi)*sind(theta)*cosd(phi))+(sind(psi)*sind(phi));
        -sind(psi)*cosd(theta) -(sind(psi)*sind(theta)*sind(phi))+(cosd(psi)*cosd(phi)) (sind(psi)*sind(theta)*cosd(phi))+(cosd(psi)*sind(phi));
        sind(theta) -cosd(theta)*sind(phi) cosd(theta)*cosd(phi)];
    D_til = obj.Dbar*D';
    obj.aTilde = [-atand(D_til(3,2)/D_til(3,3)); asind(D_til(3,1)); -atand(D_til(2,1)/D_til(1,1))];
    
        
end
