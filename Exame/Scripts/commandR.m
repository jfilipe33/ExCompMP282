function tar = commandR(t,mavsim) 
    
   if t < mavsim.t/mavsim.Ts;
       r1tar =  mavsim.ro1*cos(mavsim.omega*(t*mavsim.Ts));
       r2tar =  mavsim.ro1*sin(mavsim.omega*(t*mavsim.Ts));
       r3tar =  mavsim.ro2*(t*mavsim.Ts);
       psitar = mavsim.omega*(t*mavsim.Ts) + pi/2;
   else
       r1tar =  mavsim.ro1;
       r2tar =  0;
       r3tar =  mavsim.ro2*(mavsim.t);
       psitar = pi/2;
   end
   
   tar = [r1tar; r2tar; r3tar; psitar];
end