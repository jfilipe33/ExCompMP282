%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lista 3 - Ex. 2a)
% Author: João Filipe R. P. de A. Silva / joaofilipe@ita.br
% Affiliation: Aeronautics Institute of Technology (ITA/Brazil)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%% Simulation parameters

tf = 10 ;               % simulation time in seconds
Ts  = 0.001;            % sampling time in seconds
t = 0:Ts:tf-Ts;
ang = zeros(length(t),3);

%% MAV Properties

mav.Jb = diag([0.05 0.05 0.01]);
mav.m = 0.5;
mav.g = 9.81;


%% Attitude Controller Properties

ac.K1 = 10*eye(3);
ac.K2 = 10*eye(3);
ac.Jb = mav.Jb;

%% Variables Initialization 
ac.D = eye(3); %Attitude Matrix initialization
ac.angTar = [30 30 30]; %Euler Angles command
ac.Dbar = eye(3); %Attitude Matrix Command initialization
ac = euAng2D(ac); %Euler Angles 123 Representation
ac.DTilde = ac.Dbar*ac.D';
ac.aTilde = D2euAng(ac);
ac.alpha = [0 0 0]';
ac.alphaDot = [0 0 0]';
ac.Omega = [0 0 0]'; %Angular Velocity Matrix initialization
%ac.OmegaTilde = - ac.Omega; %Angular Velocity Error kinematics
ac.OmegaDot = [0 0 0]'; %Angular Velocity Error dynamics
states = [ac.alpha; ac.Omega];


%% Execution

for cont = 1:tf/Ts
    
    ac.OmegaTilde = -ac.Omega;
    ac.Tc = ac.Jb*ac.K1*ac.aTilde + ac.Jb*ac.K2*ac.OmegaTilde - skew(ac.Jb*ac.Omega)*ac.Omega;

    %integration using HK4
            
    k1 = Ts*dynAC(states,ac);
    k2 = Ts*dynAC(states+k1/2,ac);
    k3 = Ts*dynAC(states+k2/2,ac);         
    k4 = Ts*dynAC(states+k3,ac); 
    states  = states + k1/6 + k2/3 + k3/3 + k4/6;
    
    %updated states
            
    ac.alpha = states(1:3);
    ac.Omega = states(4:6);
    
    printAlpha(cont) = norm(ac.alpha);
    printOmega(cont) = norm(ac.Omega);
    printInput(cont) = norm(ac.Tc);
%     
    ang(cont,:) = ac.alpha;
    ac = aTildeRefresh(ac);
 
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
plot(t,printInput,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$T^c_B$','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Attitude States');
plot(t,printAlpha,'LineWidth',2.0);
plot(t,printOmega,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('Amplitude','interpreter','latex','FontSize',14);
% 
% figure; hold on; grid; box;
% title('Phase Plane');
% plot(printaTilde,printOmegaTilde,'LineWidth',2.0);
% xlabel('$$\tilde{g}$$','interpreter','latex','FontSize',14);
% ylabel('$$\dot{\tilde{g}}$$','interpreter','latex','FontSize',14);
