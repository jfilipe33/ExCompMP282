function ang = D2euAng(obj)
      
    phi = -atand(obj.D(3,2)/obj.D(3,3));
    theta = asind(obj.D(3,1));
    psi = -atand(obj.D(2,1)/obj.D(1,1));
    ang = [phi theta psi];
        
end
