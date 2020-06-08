%%%Saturation
function tosat = sat(val,LI,LS)
    tosat = val;
    for i=1:size(val,1)
        if val(i) < LI(i) 
            tosat = LI(i);
        end
        if val(i) > LS(i) 
            tosat = LS(i);
        end
    end
end
