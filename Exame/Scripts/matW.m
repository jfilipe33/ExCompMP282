function W = matW(Omega)
    W = [-skew(Omega) Omega;
        -Omega' 0];
end