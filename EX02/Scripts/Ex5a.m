%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lista 2 - Ex. 5a)
% Author: João Filipe R. P. de A. Silva / joaofilipe@ita.br
% Affiliation: Aeronautics Institute of Technology (ITA/Brazil)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%% Simulation parameters

tf = 5;               % simulation time in seconds
Ts  = 0.001;            % sampling time in seconds
t = 0:Ts:tf-Ts;
ang = zeros(length(t),3);

%% MAV Properties

mav.Jb = diag([0.05 0.05 0.01]);
mav.m = 0.5;
mav.g = 9.81;


%% Attitude Controller Properties

ac.eta = 2;
ac.lambda = -10*eye(3);
ac.Tdmax = 0.05; %Disturbance Limit
ac.C1 = -ac.lambda; %gTildeDot = -inv(C2)*C1*gTilde = Lambda*gTilde
ac.C2 = eye(3);
ac.Jb = mav.Jb;

%% Variables Initialization 
ac.D = eye(3); %Attitude Matrix initialization
ac.angTar = [30 30 30]; %Euler Angles command
ac.Dbar = eye(3); %Attitude Matrix Command initialization
ac = euAng2G(ac); %Euler Angles 123 Representation
ac.DTilde = ac.Dbar*ac.D';
ac.gTilde = (1/(1+trace(ac.DTilde))) * [ac.DTilde(2,3) - ac.DTilde(3,2); ac.DTilde(3,1) - ac.DTilde(1,3);ac.DTilde(1,2) - ac.DTilde(2,1)];
ac.Omega = [0 0 0]'; %Angular Velocity Matrix initialization
ac.OmegaTilde = - ac.Omega; %Angular Velocity Error kinematics
ac.OmegaTildeDot = [0 0 0]'; %Angular Velocity Error dynamics
%ac.gTilde = [1 1 1]'; %Attitude Error initialization (In Gibbs Vector)
ac.gTildeDot = [0 0 0]'; %Attitude Error kinematics initialization (In Gibbs Vector)
ac.M1 = eye(3);
states = [ac.gTilde; ac.OmegaTilde];


%% Execution

for cont = 1:tf/Ts
    
    ac.Td = -ac.Tdmax + 2*ac.Tdmax*rand;
    ac.M1 = ac.gTilde*ac.gTilde' + skew(ac.gTilde) + eye(3);
    ac.Ka = -ac.Tdmax - (sqrt(2)*ac.eta)/(norm(ac.C2*ac.M1*inv(ac.Jb)));
    ac.M1dot = ac.gTildeDot*ac.gTilde' + ac.gTilde*ac.gTildeDot' + skew(ac.gTildeDot);
    ac.Sa = ac.C1*ac.gTilde + ac.C2*ac.gTildeDot;
    ac.Tc = -inv(ac.C2*ac.M1*inv(ac.Jb)) * (-ac.C1*ac.M1*ac.OmegaTilde - ac.C2*ac.M1dot*ac.OmegaTilde + ac.C2*ac.M1*inv(ac.Jb)*skew(ac.Jb*ac.OmegaTilde)*ac.OmegaTilde + ac.Ka*norm(ac.C2*ac.M1*inv(ac.Jb))*sign(ac.Sa));

    %integration using HK4
            
    k1 = Ts*dynAC(states,ac);
    k2 = Ts*dynAC(states+k1/2,ac);
    k3 = Ts*dynAC(states+k2/2,ac);         
    k4 = Ts*dynAC(states+k3,ac); 
    states  = states + k1/6 + k2/3 + k3/3 + k4/6;
    
    %updated states
            
    ac.gTilde = states(1:3);
    ac.OmegaTilde = states(4:6);
    ac.gTildeDot = 0.5*(ac.gTilde*ac.gTilde' + skew(ac.gTilde) + eye(3))*ac.OmegaTilde;
    ac.OmegaTildeDot = -inv(ac.Jb)*skew(ac.Jb*ac.OmegaTilde)*ac.OmegaTilde - inv(ac.Jb)*(ac.Tc + ac.Td);
    
%     ac.Sa = ac.C1*ac.gTilde + ac.C2*ac.gTildeDot;
    printGTilde(cont) = norm(ac.gTilde);
    printGTildeDot(cont) = norm(ac.gTildeDot);
    printOmegaTilde(cont) = norm(ac.OmegaTilde);
    printInput(cont) = norm(ac.Tc);
    printS(cont) = norm(ac.Sa);
%     
    ang(cont,:) = G2euAng(ac);
 
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
title('Sliding Variable');
plot(t,printS,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$Sa$','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Attitude States');
plot(t,printGTilde,'LineWidth',2.0);
plot(t,printOmegaTilde,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('Amplitude','interpreter','latex','FontSize',14);

figure; hold on; grid; box;
title('Phase Plane');
plot(printGTilde,printGTildeDot,'LineWidth',2.0);
xlabel('$$\tilde{g}$$','interpreter','latex','FontSize',14);
ylabel('$$\dot{\tilde{g}}$$','interpreter','latex','FontSize',14);
