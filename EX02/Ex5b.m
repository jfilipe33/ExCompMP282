%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lista 2 - Ex. 5b)
% Author: João Filipe R. P. de A. Silva / joaofilipe@ita.br
% Affiliation: Aeronautics Institute of Technology (ITA/Brazil)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%% Simulation parameters

tf = 5  ;               % simulation time in seconds
Ts  = 0.001;            % sampling time in seconds
t = 0:Ts:tf-Ts;

%% MAV Properties

mav.Jb = diag([0.05 0.05 0.01]);
mav.m = 0.5;
mav.g = 9.81;
mav.e3 = [0 0 1]';

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
states = [pc.rTilde; pc.vTilde];

for cont = 1:tf/Ts
    
    pc.Fd = -pc.Fdmax + 2*pc.Fdmax*rand;
    pc.Kp = -pc.Fdmax - (mav.m*pc.eta)/(norm(pc.C4)*sqrt(2));
    pc.Sp = pc.C3*pc.rTilde + pc.C4*pc.vTilde;
    pc.Fc = mav.m*inv(pc.C4)*(pc.C3*pc.vTilde + pc.C4*mav.g*mav.e3 + pc.C4*pc.aTar - pc.Kp*(norm(pc.C4)/mav.m)*sign(pc.Sp));
    
    %integration using HK4
            
    k1 = Ts*dynPC(states,pc);
    k2 = Ts*dynPC(states+k1/2,pc);
    k3 = Ts*dynPC(states+k2/2,pc);         
    k4 = Ts*dynPC(states+k3,pc); 
    states  = states + k1/6 + k2/3 + k3/3 + k4/6;
    
    %updated states
            
    pc.rTilde = states(1:3);
    pc.vTilde = states(4:6);
    printRTilde(cont) = norm(pc.rTilde);
    printVTilde(cont) = norm(pc.vTilde);
    printInput(cont) = norm(pc.Fc);
    printS(cont) = norm(pc.Sp); 
    pos(cont,:) = pc.rTar - pc.rTilde;
 
end

figure; grid; box;
title('Position');
plot3(pos(:,1),pos(:,2),pos(:,3),'LineWidth',2.0);
xlabel('x (m)','interpreter','latex','FontSize',14);
ylabel('y (m)','interpreter','latex','FontSize',14);
zlabel('z (m)','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Input');
plot(t,printInput,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$F^c_B$','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Sliding Variable');
plot(t,printS,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$Sp$','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Position States');
plot(t,printRTilde,'LineWidth',2.0);
plot(t,printVTilde,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('Amplitude','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Phase Plane');
plot(printRTilde,printVTilde,'LineWidth',2.0);
xlabel('$$\tilde{r}$$','interpreter','latex','FontSize',14);
ylabel('$$\tilde{v}$$','interpreter','latex','FontSize',14);




