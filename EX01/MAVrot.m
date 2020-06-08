function rot = MAVrot(t,theta,tau_c)
    rot = zeros(2,1);
    rot(1) = theta(2);
    rot(2) = (1/J)*tau_c;