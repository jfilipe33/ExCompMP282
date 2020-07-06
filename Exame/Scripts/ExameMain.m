%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exame Final                                                  %
% Author: João Filipe R. P. de A. Silva / joaofilipe@ita.br    %
% Affiliation: Aeronautics Institute of Technology (ITA/Brazil)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;

%% Simulation parameters

tf = 30;               % simulation time in seconds
Ts  = 0.001;            % sampling time in seconds
t = 0:Ts:tf-Ts;
% ang = zeros(length(t),3);

%% MAV Properties

mav.Jb = diag([0.02 0.02 0.05]);
mav.m = 1;
mav.g = 9.81;
mav.e3 = [0 0 1]';
mav.delta = 22.5;
mav.l = 0.3;
mav.kf = 3.13E-5;
mav.k = 0.024;
mav.km = 1;
mav.Tm = 0.01;
mav.n = 6;
mav.fmax = 2*mav.g*mav.m/mav.n;
mav.fmin = 0.25*mav.g*mav.m/mav.n;
mav.phimax = 30;
mav.fbar = [0 0 0 0 0 0]';
mav.f = [0 0 0 0 0 0]';
mav.wbar = [0 0 0 0 0 0]';
mav.w = [0 0 0 0 0 0]';
mav.Lambda = [1 1 1 1 1 1; 
              mav.l*sind(mav.delta) -mav.l*sind(mav.delta) -mav.l -mav.l*sind(mav.delta) mav.l*sind(mav.delta) mav.l;
              -mav.l*cosd(mav.delta) -mav.l*cosd(mav.delta) 0 mav.l*cosd(mav.delta) mav.l*cosd(mav.delta) 0;
              mav.k -mav.k mav.k -mav.k mav.k -mav.k];

%% MAV Simulation Properties

mavsim.ro1 = 2;
mavsim.ro2 = 1;
mavsim.omega = 0.4*pi;
mavsim.t = 15;
mavsim.Ts = Ts;

%% Attitude Controller Properties

ac.K1 = diag([20 20 75]);
ac.K2 = diag([12 12 25]);
ac.Jb = mav.Jb;
ac.Tdmax = 0.05; %Disturbance Limit


%% Position Controller Properties

pc.K3 = diag([0.175 0.175 8]);
pc.K4 = diag([1.15 1.15 1.35]);
pc.Fdmax = 0.5; %Disturbance Limit
% pc.K3 = 1*eye(3);
% pc.K4 = 1.2*eye(3);

%% Variables Initialization 
pc.r = [2 0 0]'; %3D initial position
pc.rTar = [2 0 0]'; %3D initial position
pc.vTar = [0 0 0]'; %3D velocity command
pc.v = [0 0 0]'; %3D initial velocity
ac.a0 = [0 0 pi/2]; %Attitude initialization
ac.D = euAng2D(ac.a0); %Attitude Matrix initialization
ac.q = D2Q(ac.D); %Attitude in quaternions
ac.Dbar = eye(3); %Attitude Matrix Command initialization
ac.Omega = [0 0 0]'; %Angular Velocity Matrix initialization
states = [pc.r; pc.v; ac.q; ac.Omega; mav.w];
ac.psi = pi/2;
ac.psidot = 0;
Fc = 0;
%% Execution

for cont = 1:tf/Ts
    
    pc.rTarprev = pc.rTar;
    pc.vTarprev = pc.vTar;
    pc.tar = commandR(cont,mavsim);
    pc.rTar = pc.tar(1:3);
    rTar(cont,:) = pc.rTar;
    pc.vTar = (pc.rTar - pc.rTarprev)/mavsim.Ts;
    if norm(pc.vTar) > 10
         pc.vTar = pc.vTarprev;
    end
    pc.Fc = mav.m*(pc.K3*(pc.rTar-pc.r) + pc.K4*(pc.vTar-pc.v) + mav.g*mav.e3);
    
    ac.N = pc.Fc/norm(pc.Fc);
    ac.phi = -atan2(ac.N(2),ac.N(3));
    ac.theta = asin(ac.N(1));
    ac.psidprev = ac.psidot;
    ac.psidot = (pc.tar(4) - ac.psi)/mavsim.Ts;
    if norm(ac.psidot) > 10
         ac.psidot = ac.psidprev;
    end
    ac.psi = pc.tar(4); 
    

    ac.angTar = [ac.phi ac.theta ac.psi];
    ac.Dbar = euAng2D(ac.angTar); %Euler Angles 123 Representation
    aBarprint(cont,:) = rad2deg(D2euAng(ac.Dbar));

    ac.DTilde = ac.Dbar*ac.D';
    ac.aTilde = D2euAng(ac.DTilde);
    
    ac.OmegaTilde = [0;0;ac.psidot] - ac.Omega;

    
    ac.Td = -ac.Tdmax + 2*ac.Tdmax*rand;
    pc.Fd = -pc.Fdmax + 2*pc.Fdmax*rand;
    
    ac.Tc = ac.Jb*ac.K1*ac.aTilde + ac.Jb*ac.K2*ac.OmegaTilde - skew(ac.Jb*ac.Omega)*ac.Omega;
    
    pc.Fc = satF(pc,mav);
    ac.Tc = satT(ac,mav);
    
    mav.fbar = pinv(mav.Lambda)*[norm(pc.Fc); ac.Tc];   
    mav.wbar = sqrt(mav.fbar/mav.kf);
    mav.w = rotorTF(Ts,mav)';
    mav.f = mav.w.^2*mav.kf;
    Esf = real(mav.Lambda*mav.f);
    pc.Fc = Esf(1);
    ac.Tc = Esf(2:4);
    
    
    
%integration using HK4

k1 = Ts*dynFull(states,ac,pc,mav);
k2 = Ts*dynFull(states+k1/2,ac,pc,mav);
k3 = Ts*dynFull(states+k2/2,ac,pc,mav);         
k4 = Ts*dynFull(states+k3,ac,pc,mav); 
states  = states + k1/6 + k2/3 + k3/3 + k4/6;

%updated states

pc.r = states(1:3);
pc.v = states(4:6);
ac.q = states(7:10);
ac.q = ac.q/norm(ac.q);
ac.Omega = states(11:13);
mav.w = states(14:19);

ac.D = Q2D(ac.q);
ac.a = D2euAng(ac.D);

printR(cont) = norm(pc.r);
printV(cont) = norm(pc.v);
printInputFc(cont) = norm(pc.Fc);
% printq(cont) = norm(ac.q);
printOmega(cont) = norm(ac.Omega);
printInputTc(cont) = norm(ac.Tc);
print3(cont,:) = ac.Tc;
%     
aprint(cont,:) = rad2deg(ac.a);
pos(cont,:) = pc.r;

 end
% 
figure; hold on; grid; box;
title('Angles');
plot(t,aprint(:,1),'LineWidth',2.0);
plot(t,aprint(:,2),'LineWidth',2.0);
plot(t,aprint(:,3),'LineWidth',2.0);
plot(t,aBarprint(:,1),'LineWidth',2.0);
plot(t,aBarprint(:,2),'LineWidth',2.0);
plot(t,aBarprint(:,3),'LineWidth',2.0);
% plot(t,rad2deg(printpsi),'LineWidth',2.0);
xlabel('time (s)','interpreter','latex','FontSize',14);
ylabel('$deg$','interpreter','latex','FontSize',14);
% 
% figure; hold on; grid; box;
% title('Input');
% plot(t,printInputTc,'LineWidth',2.0);
% xlabel('time (s)','interpreter','latex','FontSize',14);
% ylabel('$T^c_B$','interpreter','latex','FontSize',14);
% 
% figure; hold on; grid; box;
% title('Torque Resultante');
% subplot(3,1,1);
% plot(t,print3(:,1),'LineWidth',2.0);
% subplot(3,1,2);
% plot(t,print3(:,2),'LineWidth',2.0);
% subplot(3,1,3);
% plot(t,print3(:,3),'LineWidth',2.0);
% xlabel('Tempo (s)','interpreter','latex','FontSize',14);
% ylabel('$T^c_B$','interpreter','latex','FontSize',14);
% 
% figure; hold on; grid; box;
% title('Attitude States');
% % plot(t,printq,'LineWidth',2.0);
% plot(t,printOmega,'LineWidth',2.0);
% xlabel('time (s)','interpreter','latex','FontSize',14);
% ylabel('Amplitude','interpreter','latex','FontSize',14);

% figure; hold on; grid; box;
% title('Posição Tridimensional X Tempo');
% subplot(3,1,1);
% plot(t,pos(:,1),'LineWidth',2.0);
% hold on;
% plot(t,rTar(:,1),'LineWidth',2.0);
% subplot(3,1,2);
% plot(t,pos(:,2),'LineWidth',2.0);
% hold on;
% plot(t,rTar(:,2),'LineWidth',2.0);
% subplot(3,1,3);
% plot(t,pos(:,3),'LineWidth',2.0);
% hold on;
% plot(t,rTar(:,3),'LineWidth',2.0);
% xlabel('Tempo (s)','interpreter','latex','FontSize',14);
% ylabel('Posição (m)','interpreter','latex','FontSize',14);
% % 
figure; grid; box;
title('Posição Tridimensional');
plot3(pos(:,1),pos(:,2),pos(:,3),'LineWidth',2.0);
xlabel('x (m)','interpreter','latex','FontSize',14);
ylabel('y (m)','interpreter','latex','FontSize',14);
zlabel('z (m)','interpreter','latex','FontSize',14);
hold on;
plot3(rTar(:,1),rTar(:,2),rTar(:,3),'LineWidth',2.0);

% figure; hold on; grid; box;
% title('3D-Position over time');
% plot(t,print2(:,1),'LineWidth',2.0);
% plot(t,print2(:,2),'LineWidth',2.0);
% plot(t,print2(:,3),'LineWidth',2.0);
% xlabel('time (s)','interpreter','latex','FontSize',14);
% ylabel('$deg$','interpreter','latex','FontSize',14);
% 
figure; hold on; grid; box;
title('Empuxo Resultante');
plot(t,printInputFc,'LineWidth',2.0);
xlabel('Tempo (s)','interpreter','latex','FontSize',14);
ylabel('$F^c_G (N)$','interpreter','latex','FontSize',14);

% figure; hold on; grid; box;
% title('Position States');
% plot(t,printR,'LineWidth',2.0);
% plot(t,printV,'LineWidth',2.0);
% xlabel('time (s)','interpreter','latex','FontSize',14);
% ylabel('Amplitude','interpreter','latex','FontSize',14);