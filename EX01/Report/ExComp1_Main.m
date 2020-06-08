%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     INSTITUTO TECNOLÓGICO DE AERONÁUTICA     %
% MP-282: Dynamic Modeling and Control of MAVs %
%          Simulation Exercise 1               %
%              2D-MAV Control                  %
%       João Filipe R. P. de A. Silva          %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Getting Ready
clc
close all
clear all

%% General Parameters

g = 9.81; %Gravity [kg*m/s²]
Ts = 0.01; %Sampling Time [s]
e1 = [1;0]; %Unitary x vector
e2 = [0;1]; %Unitary y vector
wx = [0.5 2 2 3.5 3.5]; 
wy = [2.5 2.5 1.5 1.5 2.5]; 
w = [wx;wy]'; %waypoints 

%% MAV Model

mav.m = 1; %MAV mass [kg]
mav.l = 0.2; %MAV arm length [m]
mav.v = [0;0]; %MAV speed [m/s]
mav.W = 0;  %MAV angular speed [ºx/s]
mav.J = 0.02; %MAV moment of Inertia [kg*m²]
mav.max_theta = 30; %Maximum angle of maneuver [º]
mav.r0 = [0.5;1.5]; %Initial position [m]
mav.r = mav.r0;
mav.theta0 = 0; %Initial angle [º]
mav.theta = mav.theta0;
mav.Fmax = 2*mav.m*g;
mav.Fmin = 0.75*mav.m*g;

%% Control Loop

%tspan = 0:Ts:100;
cont = 1;
Fmin = [-mav.Fmin*tand(mav.max_theta);mav.Fmin];
Fmax = [mav.Fmin*tand(mav.max_theta);mav.Fmax];
Tmax = (mav.m*g-(mav.Fmin))*mav.l;
x = [mav.r;mav.v;mav.theta;mav.W];
K1 = 6;
K2 = 3;
K3 = eye(2)*0.2;
K4 = eye(2)*0.8;
Kg = 1.5*eye(2);

waycount = 1;
while waycount < 6
   
   mav.r_= mav.r + Kg*((w(waycount,:)'-mav.r));
   %mav.r_= mav.r + Ts*0.5*((w(waycount,:)'-mav.r)/norm(w(waycount,:)'-mav.r));
   mav.D = [cosd(mav.theta) -sind(mav.theta);sind(mav.theta) cosd(mav.theta)]';
   mav.Fc = mav.m*K3*(mav.r_-mav.r)-mav.m*K4*mav.v+mav.m*g*e2;
   mav.Fc = sat(mav.Fc,Fmin,Fmax);
   mav.theta_ = atand(e1'*mav.Fc)/(e2'*mav.Fc);
   mav.Tc = (mav.J*K1*(mav.theta_ - mav.theta) - mav.J*K2*mav.W);
   mav.Tc = sat(mav.Tc,-Tmax,Tmax);
   mav.f = CA(mav.Fc,mav.Tc,mav.l);

%integration using HK4
            
   k1 = Ts*dyn(x,mav);
   k2 = Ts*dyn(x+k1/2,mav);
   k3 = Ts*dyn(x+k2/2,mav);            
   k4 = Ts*dyn(x+k3,mav); 
   x  = x + k1/6 + k2/3 + k3/3 + k4/6;

%updated states
            
   mav.r = x(1:2);
   mav.v = x(3:4);
   mav.theta = x(5);
   mav.W = x(6);
   
   rp(cont,1) = mav.r(1);
   rp(cont,2) = mav.r(2);
   fp(cont,1) = mav.f(1);
   fp(cont,2) = mav.f(2);
   tp(cont) = mav.theta;
   tbp(cont) = mav.theta_;
   Tp(cont) = mav.Tc;
   Fp(cont,1) = mav.Fc(1);
   Fp(cont,2) = mav.Fc(2);
   wp(cont) = mav.W;
   cont = cont + 1;
   
   if bounded(mav.r(1),mav.r(2),w(waycount,1),w(waycount,2)) == 1;
       waycount = waycount + 1;
   end
   
   mav.r
end
stencil()
time = 0:Ts:Ts*(cont-2);
plot(rp(:,1),rp(:,2));
figure 
plot(time,fp(:,1),time,fp(:,2));
figure
plot(time,tp,time,tbp);
figure
plot(time,wp)
% figure
% plot(1:cont-1,Fp)