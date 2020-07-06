function Tout = satT(obja, objm)
   
   Tout = obja.Tc;
   zetamin = objm.fmin;
   zetamax = objm.fmax - objm.fmin;
   Txmax = objm.l*(zetamax - zetamin)*(1+2*sind(objm.delta));
   Txmin = -Txmax;
   Tymax = 2*objm.l*(zetamax - zetamin)*(cosd(objm.delta));
   Tymin = -Tymax;
   Tzmax = 3*objm.k*(zetamax - zetamin);
   Tzmin = -Tzmax;
   if obja.Tc(1)<Txmin, Tout(1) = Txmin; end
   if obja.Tc(1)>Txmax, Tout(1) = Txmax; end
   if obja.Tc(2)<Tymin, Tout(2) = Tymin; end
   if obja.Tc(2)>Tymax, Tout(2) = Tymax; end
   if obja.Tc(3)<Tzmin, Tout(3) = Tzmin; end
   if obja.Tc(3)>Tzmin, Tout(3) = Tzmax; end
   
end