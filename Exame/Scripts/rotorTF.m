function w = rotorTF(Ts,obj)
   
    for i=1:obj.n
                
        w(i) = exp(-Ts/obj.Tm)*(obj.w(i)-obj.km*obj.wbar(i))+ obj.km*obj.wbar(i);
                     
    end
   
end