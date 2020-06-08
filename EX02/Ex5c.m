%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lista 2 - Ex. 5c)
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
%ang = zeros(length(t),3);
sigmo = 0.02;

%% MAV Properties

mav.Jb = diag([0.05 0.05 0.01]);
mav.m = 0.5;
mav.g = 9.81;
mav.e3 = [0 0 1]';

%% Attitude Controller Properties

ac.eta = 20;
ac.lambda = -10*eye(3);
ac.Tdmax = 0.05; %Disturbance Limit
ac.C1 = -ac.lambda; %gTildeDot = -inv(C2)*C1*gTilde = Lambda*gTilde
ac.C2 = eye(3);
ac.Jb = mav.Jb;

%% Position Controller Properties

pc.eta = 2;
pc.lambda = -1*eye(3);
pc.Fdmax = 0.5; %Disturbance Limit
pc.C3 = -pc.lambda; %rTildeDot = -inv(C4)*C3*rTilde = Lambda*rTilde
pc.C4 = eye(3);
pc.Jb = mav.Jb;

%% Variables Initialization 

pc.rTar = [1 1 1]'; %3D position command
pc.r = [0 0 0]'; %3D initial position
pc.rTilde = pc.rTar - pc.r; %3D position error
pc.vTar = [0 0 0]'; %3D velocity command
pc.v = [0 0 0]'; %3D initial velocity
pc.vTilde = pc.vTar - pc.v; %3D velocity error
pc.aTar = [0 0 0]'; %3D acceleration command

ac.D = eye(3); %Attitude Matrix initialization
% ac.angTar = [30 30 30]; %Euler Angles command
ac.Dbar = eye(3); %Attitude Matrix Command initialization
%ac.DTilde = ac.Dbar*ac.D';
%ac.gTilde = (1/(1+trace(ac.DTilde))) * [ac.DTilde(2,3) - ac.DTilde(3,2); ac.DTilde(3,1) - ac.DTilde(1,3);ac.DTilde(1,2) - ac.DTilde(2,1)];
ac.Omega = [0 0 0]'; %Angular Velocity Matrix initialization
ac.OmegaTilde = - ac.Omega; %Angular Velocity Error kinematics
ac.OmegaTildeDot = [0 0 0]'; %Angular Velocity Error dynamics
%ac.gTilde = [1 1 1]'; %Attitude Error initialization (In Gibbs Vector)
ac.gTildeDot = [0 0 0]'; %Attitude Error kinematics initialization (In Gibbs Vector)
ac.M1 = eye(3);


%% Execution

for cont = 1:tf/Ts
    
    pc.Fd = -pc.Fdmax + 2*pc.Fdmax*rand;
    pc.Kp = -pc.Fdmax - (mav.m*pc.eta)/(norm(pc.C4)*sqrt(2));
    pc.Sp = pc.C3*pc.rTilde + pc.C4*pc.vTilde;
    %Comentar linha abaixo para usar sigmoide
    %pc.Fc = mav.m*inv(pc.C4)*(pc.C3*pc.vTilde + pc.C4*mav.g*mav.e3 + pc.C4*pc.aTar - pc.Kp*(norm(pc.C4)/mav.m)*(sign(pc.Sp)));
    %Comentar linha abaixo para usar sign(s)
    pc.Fc = mav.m*inv(pc.C4)*(pc.C3*pc.vTilde + pc.C4*mav.g*mav.e3 + pc.C4*pc.aTar - pc.Kp*(norm(pc.C4)/mav.m)*(pc.Sp/(norm(pc.Sp) + sigmo)));
    %pc = fSat(mav,pc);
    
    ac.N = pc.Fc/norm(pc.Fc);
    ac.phi = -atand(ac.N(2)/ac.N(3));
    ac.theta = asind(ac.N(1));
    ac.psi = 0;
    ac.angTar = [ac.phi ac.theta ac.psi];
    ac = euAng2G(ac); %Euler Angles 123 Representation
    ac.DTilde = ac.Dbar*ac.D';
    ac.gTilde = (1/(1+trace(ac.DTilde))) * [ac.DTilde(2,3) - ac.DTilde(3,2); ac.DTilde(3,1) - ac.DTilde(1,3);ac.DTilde(1,2) - ac.DTilde(2,1)];
    states = [pc.rTilde; pc.vTilde; ac.gTilde; ac.OmegaTilde];
    
    ac.Td = -ac.Tdmax + 2*ac.Tdmax*rand;
    ac.M1 = ac.gTilde*ac.gTilde' + skew(ac.gTilde) + eye(3);
    ac.Ka = -ac.Tdmax - (sqrt(2)*ac.eta)/(norm(ac.C2*ac.M1*inv(ac.Jb)));
    ac.M1dot = ac.gTildeDot*ac.gTilde' + ac.gTilde*ac.gTildeDot' + skew(ac.gTildeDot);
    ac.Sa = ac.C1*ac.gTilde + ac.C2*ac.gTildeDot;
    %Comentar linha abaixo para usar sigmoide
    %ac.Tc = -inv(ac.C2*ac.M1*inv(ac.Jb)) * (-ac.C1*ac.M1*ac.OmegaTilde - ac.C2*ac.M1dot*ac.OmegaTilde + ac.C2*ac.M1*inv(ac.Jb)*skew(ac.Jb*ac.OmegaTilde)*ac.OmegaTilde + ac.Ka*norm(ac.C2*ac.M1*inv(ac.Jb))*sign(ac.Sa));
    %Comentar linha abaixo para usar sign(s)
    ac.Tc = -inv(ac.C2*ac.M1*inv(ac.Jb)) * (-ac.C1*ac.M1*ac.OmegaTilde - ac.C2*ac.M1dot*ac.OmegaTilde + ac.C2*ac.M1*inv(ac.Jb)*skew(ac.Jb*ac.OmegaTilde)*ac.OmegaTilde + ac.Ka*norm(ac.C2*ac.M1*inv(ac.Jb))*(ac.Sa/(norm(ac.Sa) + sigmo)));

    %integration using HK4
            
    k1 = Ts*dynFull(states,ac,pc);
    k2 = Ts*dynFull(states+k1/2,ac,pc);
    k3 = Ts*dynFull(states+k2/2,ac,pc);         
    k4 = Ts*dynFull(states+k3,ac,pc); 
    states  = states + k1/6 + k2/3 + k3/3 + k4/6;
    
    %updated states
            
    pc.rTilde = states(1:3);
    pc.vTilde = states(4:6);
    ac.gTilde = states(7:9);
    ac.OmegaTilde = states(10:12);
    
    ac.gTildeDot = 0.5*(ac.gTilde*ac.gTilde' + skew(ac.gTilde) + eye(3))*ac.OmegaTilde;
    ac.OmegaTildeDot = -inv(ac.Jb)*skew(ac.Jb*ac.OmegaTilde)*ac.OmegaTilde - inv(ac.Jb)*(ac.Tc + ac.Td);
    
    ac = G2D(ac);
    
%     ac.Sa = ac.C1*ac.gTilde + ac.C2*ac.gTildeDot;
    printRTilde(cont) = norm(pc.rTilde);
    printVTilde(cont) = norm(pc.vTilde);
    printGTilde(cont) = norm(ac.gTilde);
    printGTildeDot(cont) = norm(ac.gTildeDot);
    printOmegaTilde(cont) = norm(ac.OmegaTilde);
    printInputT(cont) = norm(ac.Tc);
    printInputF(cont) = norm(pc.Fc);
    printSa(cont) = norm(ac.Sa);
    printSp(cont) = norm(pc.Sp);
    pos(cont,:) = pc.rTar - pc.rTilde;
    ang(cont,:) = D2euAng(ac);
%     
    
 
end

figure; grid; box;
title('Position');
plot3(pos(:,1),pos(:,2),pos(:,3),'LineWidth',2.0);
xlabel('x (m)','interpreter','latex','FontSize',14);
ylabel('y (m)','interpreter','latex','FontSize',14);
zlabel('z (m)','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Thrust Command');
plot(t,printInputF,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$F^c_B$','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Position Sliding Variable');
plot(t,printSp,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$Sp$','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Position States');
plot(t,printRTilde,'LineWidth',2.0);
plot(t,printVTilde,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('Amplitude','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Position Phase Plane');
plot(printRTilde,printVTilde,'LineWidth',2.0);
xlabel('$$\tilde{r}$$','interpreter','latex','FontSize',14);
ylabel('$$\tilde{v}$$','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Angles');
plot(t,ang(:,1),'LineWidth',2.0);
plot(t,ang(:,2),'LineWidth',2.0);
plot(t,ang(:,3),'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$deg$','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Torque Command');
plot(t,printInputT,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$T^c_B$','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Attitude Sliding Variable');
plot(t,printSa,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$Sa$','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Attitude States');
plot(t,printGTilde,'LineWidth',2.0);
plot(t,printOmegaTilde,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('Amplitude','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Attitude Phase Plane');
plot(printGTilde,printGTildeDot,'LineWidth',2.0);
xlabel('$$\tilde{g}$$','interpreter','latex','FontSize',14);
ylabel('$$\dot{\tilde{g}}$$','interpreter','latex','FontSize',14);
