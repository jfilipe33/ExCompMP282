function obj = fSat(mav,obj)

    fcmin = 0.5*mav.m*mav.g;
    fcmax = 2*mav.m*mav.g;
    
    if(obj.Fc(3) > fcmax)
        obj.Fc(3) = fcmax;
    elseif(obj.Fc(3) < fcmin)
        obj.Fc(3) = fcmin;
    else
        obj.Fc(3) = obj.Fc(3);
    end
    
end