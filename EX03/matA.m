function A = matA(ang)
    A = [cosd(ang(3))/cosd(ang(2)) -sind(ang(3))/cosd(ang(2)) 0;
        sind(ang(3)) cosd(ang(3)) 0;
        -(cosd(ang(3))*sind(ang(2)))/cosd(ang(2)) (sind(ang(3))*sind(ang(2)))/cosd(ang(2)) 1];
end