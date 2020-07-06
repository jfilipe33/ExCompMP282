function Fout = satF(objp, objm)

   fxmax = objm.m*objm.g*tand(objm.phimax);
   fymax = fxmax;
   fxmin = -fxmax;
   fymin = -fymax;
   fzmin = objm.fmin*objm.n;
   fzmax = objm.fmax*objm.n;
   Fout = objp.Fc;
   
   if objp.Fc(1)<fxmin, Fout(1) = fxmin; end
   
   if objp.Fc(1)>fxmax, Fout(1) = fxmax; end
   
   if objp.Fc(2)<fymin, Fout(2) = fymin; end
   
   if objp.Fc(2)>fymax, Fout(2) = fymax; end
   
   if objp.Fc(3)<fzmin, Fout(3) = fzmin; end
   
   if objp.Fc(3)>fzmax, Fout(3) = fzmax; end
   
end