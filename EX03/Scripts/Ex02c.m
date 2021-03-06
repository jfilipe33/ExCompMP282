%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lista 3 - Ex. 2c)
% Author: João Filipe R. P. de A. Silva / joaofilipe@ita.br
% Affiliation: Aeronautics Institute of Technology (ITA/Brazil)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%% Simulation parameters

tf = 10;               % simulation time in seconds
Ts  = 0.001;            % sampling time in seconds
t = 0:Ts:tf-Ts;
ang = zeros(length(t),3);

%% MAV Properties

mav.Jb = diag([0.05 0.05 0.01]);
mav.m = 0.5;
mav.g = 9.81;
mav.e3 = [0 0 1]';


%% Attitude Controller Properties

ac.K1 = 10*eye(3);
ac.K2 = 10*eye(3);
ac.Jb = mav.Jb;


%% Position Controller Properties

pc.K3 = 1.2*eye(3);
pc.K4 = 2*eye(3);

%% Variables Initialization 
pc.rTar = [1 1 1]'; %3D position command
pc.r = [0 0 0]'; %3D initial position
pc.vTar = [0 0 0]'; %3D velocity command
pc.v = [0 0 0]'; %3D initial velocity
pc.aTar = [0 0 0]'; %3D acceleration command
ac.D = eye(3); %Attitude Matrix initialization
% ac.angTar = [30 30 30]; %Euler Angles command
ac.Dbar = eye(3); %Attitude Matrix Command initialization
ac.Omega = [0 0 0]'; %Angular Velocity Matrix initialization
ac.alpha = [0 0 0]';
ac.alphaDot = [0 0 0]';
states = [pc.r; pc.v; ac.alpha; ac.Omega];

%% Execution

for cont = 1:tf/Ts
    
    pc.Fc = mav.m*(pc.K3*(pc.rTar-pc.r) + pc.K4*(pc.vTar-pc.v) + mav.g*mav.e3);
    ac.N = pc.Fc/norm(pc.Fc);
    ac.phi = -atand(ac.N(2)/ac.N(3));
    ac.theta = asind(ac.N(1));
    ac.psi = 30;
    ac.angTar = [ac.phi ac.theta ac.psi];
    ac = euAng2D(ac); %Euler Angles 123 Representation
    ac.DTilde = ac.Dbar*ac.D';
    ac.aTilde = D2euAng(ac);
    ac.OmegaTilde = -ac.Omega;
    
    ac.Tc = ac.Jb*ac.K1*ac.aTilde + ac.Jb*ac.K2*ac.OmegaTilde - skew(ac.Jb*ac.Omega)*ac.Omega;
    
    %integration using HK4
            
    k1 = Ts*dynFull(states,ac,pc);
    k2 = Ts*dynFull(states+k1/2,ac,pc);
    k3 = Ts*dynFull(states+k2/2,ac,pc);         
    k4 = Ts*dynFull(states+k3,ac,pc); 
    states  = states + k1/6 + k2/3 + k3/3 + k4/6;
    
    %updated states
    
    pc.r = states(1:3);
    pc.v = states(4:6);

    ac.alpha = states(7:9);
    ac.Omega = states(10:12);

    ac = DRefresh(ac);
    
    printR(cont) = norm(pc.r);
    printV(cont) = norm(pc.v);
    printInputFc(cont) = norm(pc.Fc);
    printAlpha(cont) = norm(ac.alpha);
    printOmega(cont) = norm(ac.Omega);
    printInputTc(cont) = norm(ac.Tc);
%     
    ang(cont,:) = ac.alpha;
    pos(cont,:) = pc.r;
 
end

figure; hold on; grid; box;
title('Angles');
plot(t,ang(:,1),'LineWidth',2.0);
plot(t,ang(:,2),'LineWidth',2.0);
plot(t,ang(:,3),'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$deg$','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Input');
plot(t,printInputTc,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$T^c_B$','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Attitude States');
plot(t,printAlpha,'LineWidth',2.0);
plot(t,printOmega,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('Amplitude','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('3D-Position over time');
plot(t,pos(:,1),'LineWidth',2.0);
plot(t,pos(:,2),'LineWidth',2.0);
plot(t,pos(:,3),'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$deg$','interpreter','latex','FontSize',14);

figure; grid; box;
title('Position');
plot3(pos(:,1),pos(:,2),pos(:,3),'LineWidth',2.0);
xlabel('x (m)','interpreter','latex','FontSize',14);
ylabel('y (m)','interpreter','latex','FontSize',14);
zlabel('z (m)','interpreter','latex','FontSize',14);


figure; hold on; grid; box;
title('Input');
plot(t,printInputFc,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$F^c_B$','interpreter','latex','FontSize',14);


figure; hold on; grid; box;
title('Position States');
plot(t,printR,'LineWidth',2.0);
plot(t,printV,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('Amplitude','interpreter','latex','FontSize',14);