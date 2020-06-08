function ang = D2euAng(obj)
      
    phi = -atand(obj.DTilde(3,2)/obj.DTilde(3,3));
    theta = asind(obj.DTilde(3,1));
    psi = -atand(obj.DTilde(2,1)/obj.DTilde(1,1));
    ang = [phi theta psi]';
        
end
