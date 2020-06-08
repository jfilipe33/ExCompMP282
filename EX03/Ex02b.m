%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lista 3 - Ex. 2b)
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

%% MAV Properties

mav.Jb = diag([0.05 0.05 0.01]);
mav.m = 0.5;
mav.g = 9.81;
mav.e3 = [0 0 1]';


%% Position Controller Properties

pc.K3 = 25*eye(3);
pc.K4 = 10*eye(3);
pc.Jb = mav.Jb;

%% Variables Initialization 
pc.rTar = [1 1 1]'; %3D position command
pc.r = [0 0 0]'; %3D initial position
pc.v = [0 0 0]'; %3D initial velocity
pc.vTar = [0 0 0]'; %3D velocity command
states = [pc.r; pc.v];


%% Execution

for cont = 1:tf/Ts
    
    pc.Fc = mav.m*(pc.K3*(pc.rTar-pc.r) + pc.K4*(pc.vTar-pc.v) + mav.g*mav.e3);

    %integration using HK4
            
    k1 = Ts*dynPC(states,pc);
    k2 = Ts*dynPC(states+k1/2,pc);
    k3 = Ts*dynPC(states+k2/2,pc);         
    k4 = Ts*dynPC(states+k3,pc); 
    states  = states + k1/6 + k2/3 + k3/3 + k4/6;
    
    %updated states
            
    pc.r = states(1:3);
    pc.v = states(4:6);
    pc.rDot = pc.v;
    pc.vDot = (1/mav.m)*pc.Fc - mav.g*mav.e3;
    
%     ac.Sa = ac.C1*ac.gTilde + ac.C2*ac.gTildeDot;
    printR(cont) = norm(pc.r);
    printV(cont) = norm(pc.v);
    printInput(cont) = norm(pc.Fc);
%     
    pos(cont,:) = pc.r;
%   pc = aTildeRefresh(pc);
 
end

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
plot(t,printInput,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$F^c_B$','interpreter','latex','FontSize',14);


figure; hold on; grid; box;
title('Position States');
plot(t,printR,'LineWidth',2.0);
plot(t,printV,'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('Amplitude','interpreter','latex','FontSize',14);