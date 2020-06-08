function states = dynAC(st,obj) 
    
    g_Til = st(1:3);
    W_Til = st(4:6);

    %obj.gTildeDot = 0.5*obj.M1*obj.OmegaTilde;
    %obj.OmegaTildeDot = -inv(obj.Jb)*skew(obj.Jb*obj.OmegaTilde)*obj.OmegaTilde - inv(obj.Jb)*(obj.Tc + obj.Td);
    
    states = [0.5*(g_Til*g_Til' + skew(g_Til) + eye(3))*W_Til;
              -inv(obj.Jb)*skew(obj.Jb*W_Til)*W_Til - inv(obj.Jb)*(obj.Tc + obj.Td)];
          
%     obj.gTildeDot = states(1:3);
%     obj.OmegaTildeDot = states(4:6);

end