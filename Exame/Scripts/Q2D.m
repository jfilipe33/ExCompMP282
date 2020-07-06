function D = Q2D(q)
   
    eps = q(1:3);
    eta = q(4); 

    D  = (eta^2-eps'*eps)*eye(3) + 2*eps*(eps') - 2*eta*skew(eps);
end