clc; clear; close all

set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLineWidth',1.2)
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesGridAlpha',0.3)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
lineStyles   = {'--','-',':','-.'};
markerStyles = {'.','o','x'};

%% Simulation Setup

% Define orbit elements
a = 10000; % km
e = 0.001;
i = deg2rad(40); % deg -> rad
Omega = deg2rad(80); % deg -> rad
w = deg2rad(40); % deg -> rad
nu0 = 0; % rad

% Define Earth Grav Model
constants.mu = 3.986004415e5; % km^3/s^2
constants.J2 = 1.0826269e-3;
constants.J3 = -2.5324e-6;
constants.RE = 6378; % km
constants.wE = (2*pi)/24 * 1/3600; % rad/s
constants.theta0 = deg2rad(122);

% Orbital period
n = sqrt(constants.mu/a^3);
T = (2*pi)/n;

% Time vector
tspan = [0 15*T];

% Station Parameters
s1_lla = [-35.398333, 148.981944, 0]; % deg, deg, m
s2_lla = [40.427222, 355.739444, 0]; % deg, deg, m
s3_lla = [35.247164, 243.205, 0]; % deg, deg, m
stations.lla = [s1_lla; s2_lla; s3_lla];
stations.station_ids = 1:length(stations.lla);
stations.el_mask = deg2rad(10); % deg -> rad

% Convert lla to ECEF
stations.Rs = lla2ecef(stations.lla)./1000;

% Get initial conditions
[r0,v0] = oe2rv(constants.mu,a,e,Omega,i,w,nu0);
perturbation = [1 1 1 1e-3 1e-3 1e-3]';
X0 = [r0;v0] + perturbation;

% ode45 calls for truth data
options = odeset('RelTol',1e-11,'AbsTol',1e-11);
[tJ2,XnomJ2] = ode45(@(t,x) orbitEOM_J2(t,x,constants),tspan,X0,options);
[t,Xnom] = ode45(@(t,x) orbitEOM_J2_J3(t,x,constants),tspan,X0,options);

% Generate measurements
measurements = genMeasurements(t,stations,constants,Xnom);
measurementsJ2 = genMeasurements(tJ2,stations,constants,XnomJ2);
t_meas = measurements(:,1);
t_measJ2 = measurementsJ2(:,1);
station_id = measurements(:,2);

% Measurement covariance matrix
sigma_r = 1e-3; % [km]
sigma_rdot = 1e-6; % [km/s] 
R = diag([sigma_r,sigma_rdot].^2);

% Add noise to measurements
noisy_measurements = zeros(size(measurements));
noisy_measurementsJ2 = zeros(size(measurementsJ2));
for i = 1:length(t_meas)
    M = measurements(i,:);
    MJ2 = measurementsJ2(i,:);
    n = size(M,1);
    M(3) = M(3) + sigma_r*randn(n,1);
    M(4) = M(4) + sigma_rdot*randn(n,1);
    MJ2(3) = MJ2(3) + sigma_r*randn(n,1);
    MJ2(4) = MJ2(4) + sigma_rdot*randn(n,1);
    noisy_measurementsJ2(i,:) = MJ2;
end

% Redo ode with only measurement time steps
[tJ2,XnomJ2] = ode45(@(t,x) orbitEOM_J2(t,x,constants),t_measJ2,X0,options);
[t,Xnom] = ode45(@(t,x) orbitEOM_J2_J3(t,x,constants),t_meas,X0,options);

%% SNC Implementation

% Initial priori covariance and state error
Pbar0 = diag([1,1,1,1e-3,1e-3,1e-3].^2);
xbar0 = [1 1 1 1e-3 1e-3 1e-3]';

% Save data
save('simulation_dataJ2J3_test.mat','t','Xnom',"noisy_measurements","R",'Pbar0',"xbar0","constants","stations")
save('simulation_dataJ2_test.mat','tJ2','XnomJ2',"noisy_measurementsJ2","R",'Pbar0',"xbar0","constants","stations")

% Iterate for different Q's
% sigma_xyz = 10.^(-15:-6);
% for i = 1:length(sigma_xyz)
% 
%     % Process noise covariance
%     Q = diag([sigma_xyz(i), sigma_xyz(i), sigma_xyz(i)].^2);
% 
%     % LKF
%     num_iterations = 1;
%     [XhatLKF,PLKF,yLKF,yhatLKF] = orbitLKF(t,xbar0,Pbar0,Q,"ECI",noisy_measurements,R,Xnom,constants,stations,num_iterations);
%     % plot_filter_diagnostics(t,Xnom,XhatLKF(:,:,end),PLKF(:,:,:,end),yLKF(:,:,end),yhatLKF(:,:,end),R,station_id,"LKF","")
% 
%     % Post-fit RMS
%     rms_pfLKF(:,i) = rms(yhatLKF);
% 
%     % State errors
%     e = XhatLKF - Xnom;
%     er = e(:,1:3); 
%     ev = e(:,4:6); 
%     er3 = sqrt(sum(er.^2,2));
%     ev3 = sqrt(sum(ev.^2,2));
% 
%     rms_rLKF(i) = rms(er3);
%     rms_vLKF(i) = rms(ev3);
% 
%     % EKF
%     LKFinit = 100;
%     [XhatEKF,PEKF,yEKF,yhatEKF] = orbitEKF(t,xbar0,Pbar0,Q,"ECI",noisy_measurements,R,Xnom,constants,stations,LKFinit);
%     % plot_filter_diagnostics(t,Xnom,XhatEKF,PEKF,yEKF,yhatEKF,R,station_id,"EKF","")
% 
%     % Post-fit RMS
%     rms_pfEKF(:,i) = rms(yhatEKF);
% 
%     % State errors
%     e = XhatEKF - Xnom;
%     er = e(:,1:3); 
%     ev = e(:,4:6); 
%     er3 = sqrt(sum(er.^2,2));
%     ev3 = sqrt(sum(ev.^2,2));
% 
%     rms_rEKF(i) = rms(er3);
%     rms_vEKF(i) = rms(ev3);
% end
% 
% % RMS vs sigma plots
% figure(); % LKF
% subplot(2,1,1)
% plot(sigma_xyz,rms_rLKF)
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('3D Position RMS [km]')
% subplot(2,1,2)
% plot(sigma_xyz,rms_vLKF)
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('3D Velocity RMS [km/s]')
% sgtitle('3D State Error RMS Values vs \sigma{x,y,z} for a LKF')
% 
% figure(); % LKF
% subplot(2,1,1)
% plot(sigma_xyz,rms_pfLKF(1,:))
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('Range Post-fit Residual RMS [km]')
% subplot(2,1,2)
% plot(sigma_xyz,rms_pfLKF(2,:))
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('Range Rate Post-fit Residual RMS [km/s]')
% sgtitle('Post-fit Residual Error RMS Values vs \sigma{x,y,z} for a LKF')
% 
% figure(); % EKF
% subplot(2,1,1)
% plot(sigma_xyz,rms_rEKF)
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('3D Position RMS [km]')
% subplot(2,1,2)
% plot(sigma_xyz,rms_vEKF)
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('3D Velocity RMS [km/s]')
% sgtitle('3D State Error RMS Values vs \sigma{x,y,z} for an EKF')
% 
% figure(); % EKF
% subplot(2,1,1)
% plot(sigma_xyz,rms_pfEKF(1,:))
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('Range Post-fit Residual RMS [km]')
% subplot(2,1,2)
% plot(sigma_xyz,rms_pfEKF(2,:))
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('Range Rate Post-fit Residual RMS [km/s]')
% sgtitle('Post-fit Residual Error RMS Values vs \sigma{x,y,z} for an EKF')

% Best Q
sigma_xyz = 1e-7;
Q = diag([sigma_xyz, sigma_xyz, sigma_xyz].^2);

% LKF
num_iterations = 1;
% [XhatLKF,PLKF,yLKF,yhatLKF] = orbitLKF(t,xbar0,Pbar0,Q,"ECI",noisy_measurements,R,Xnom,constants,stations,num_iterations);
% plot_filter_diagnostics(t,Xnom,XhatLKF(:,:,end),PLKF(:,:,:,end),yLKF(:,:,end),yhatLKF(:,:,end),R,station_id,"LKF","")

% EKF
LKFinit = 100;
[XhatEKF,PEKF,yEKF,yhatEKF] = orbitEKF(t,xbar0,Pbar0,Q,"ECI",noisy_measurements,R,Xnom,constants,stations,LKFinit);
plot_filter_diagnostics(t,Xnom,XhatEKF,PEKF,yEKF,yhatEKF,R,station_id,"EKF","")

% Define in RIC frame
% Q = diag([1e-10 1e-9 1e-7].^2);
% [XhatLKF,PLKF,yLKF,yhatLKF] = orbitLKF(t,xbar0,Pbar0,Q,"RIC",noisy_measurements,R,Xnom,constants,stations,num_iterations);
% plot_filter_diagnostics(t,Xnom,XhatLKF(:,:,end),PLKF(:,:,:,end),yLKF(:,:,end),yhatLKF(:,:,end),R,station_id,"LKF","")
% 
% [XhatEKF,PEKF,yEKF,yhatEKF] = orbitEKF(t,xbar0,Pbar0,Q,"RIC",noisy_measurements,R,Xnom,constants,stations,LKFinit);
% plot_filter_diagnostics(t,Xnom,XhatEKF,PEKF,yEKF,yhatEKF,R,station_id,"EKF","")

%% DMC Implementation

% Truth J3 accelerations
for i = 1:length(t)    
    x = Xnom(i,1);
    y = Xnom(i,2);
    z = Xnom(i,3);
    r = sqrt(x^2 + y^2 + z^2);
    uJ3_coeff = (0.5*constants.J3*constants.mu*constants.RE^3)/r^5;
    aJ3x = uJ3_coeff*(5*x/r*(7*z^3/r^3 - 3*z/r));
    aJ3y = uJ3_coeff*(5*y/r*(7*z^3/r^3 - 3*z/r));
    aJ3z = uJ3_coeff*(35*z^4/r^4 - 30*z^2/r^2+3);
    aJ3(i,:) = [aJ3x aJ3y aJ3z];
end
XnomDMC = [Xnom, aJ3];

% DMC Tuning
constants.tau = T/30;

% sigma_xyz_DMC = 10.^(-15:-5);
% for i = 1:length(sigma_xyz_DMC)
%     Q_DMC = diag([sigma_xyz_DMC(i), sigma_xyz_DMC(i), sigma_xyz_DMC(i)].^2);
% 
%     % LKF
%     [XhatLKF_DMC,PLKF_DMC,yLKF_DMC,yhatLKF_DMC] = orbitLKF_DMC(t,xbar0DMC,Pbar0DMC,Q_DMC,noisy_measurements,R,XnomDMC,constants,stations,num_iterations);
%     % plot_filter_diagnostics(t,XnomDMC,XhatLKF_DMC(:,:,end),PLKF_DMC(:,:,:,end),yLKF_DMC(:,:,end),yhatLKF_DMC(:,:,end),R,station_id,"LKF with DMC","")
% 
%     % Post-fit RMS
%     rms_pfLKF(:,i) = rms(yhatLKF_DMC);
% 
%     % State errors
%     e = XhatLKF_DMC - XnomDMC;
%     er = e(:,1:3); 
%     ev = e(:,4:6); 
%     er3 = sqrt(sum(er.^2,2));
%     ev3 = sqrt(sum(ev.^2,2));
% 
%     rms_rLKF(i) = rms(er3);
%     rms_vLKF(i) = rms(ev3);
% 
%     % EKF
%     [XhatEKF_DMC,PEKF_DMC,yEKF_DMC,yhatEKF_DMC] = orbitEKF_DMC(t,xbar0DMC,Pbar0DMC,Q_DMC,noisy_measurements,R,XnomDMC,constants,stations,LKFinit);
%     % plot_filter_diagnostics(t,XnomDMC,XhatEKF_DMC(:,:,end),PEKF_DMC(:,:,:,end),yEKF_DMC(:,:,end),yhatEKF_DMC(:,:,end),R,station_id,"EKF with DMC","")
% 
%     % Post-fit RMS
%     rms_pfEKF(:,i) = rms(yhatEKF_DMC);
% 
%     % State errors
%     e = XhatEKF_DMC - XnomDMC;
%     er = e(:,1:3); 
%     ev = e(:,4:6); 
%     er3 = sqrt(sum(er.^2,2));
%     ev3 = sqrt(sum(ev.^2,2));
% 
%     rms_rEKF(i) = rms(er3);
%     rms_vEKF(i) = rms(ev3);
% 
% end
% 
% % RMS vs sigma plots
% figure(); % LKF
% subplot(2,1,1)
% plot(sigma_xyz_DMC,rms_rLKF)
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('3D Position RMS [km]')
% subplot(2,1,2)
% plot(sigma_xyz_DMC,rms_vLKF)
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('3D Velocity RMS [km/s]')
% sgtitle('3D State Error RMS Values vs \sigma{x,y,z} for a LKF')
% 
% figure(); % LKF
% subplot(2,1,1)
% plot(sigma_xyz_DMC,rms_pfLKF(1,:))
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('Range Post-fit Residual RMS [km]')
% subplot(2,1,2)
% plot(sigma_xyz_DMC,rms_pfLKF(2,:))
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('Range Rate Post-fit Residual RMS [km/s]')
% sgtitle('Post-fit Residual Error RMS Values vs \sigma{x,y,z} for a LKF')
% 
% figure(); % EKF
% subplot(2,1,1)
% plot(sigma_xyz_DMC,rms_rEKF)
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('3D Position RMS [km]')
% subplot(2,1,2)
% plot(sigma_xyz_DMC,rms_vEKF)
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('3D Velocity RMS [km/s]')
% sgtitle('3D State Error RMS Values vs \sigma{x,y,z} for an EKF')
% 
% figure(); % EKF
% subplot(2,1,1)
% plot(sigma_xyz_DMC,rms_pfEKF(1,:))
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('Range Post-fit Residual RMS [km]')
% subplot(2,1,2)
% plot(sigma_xyz_DMC,rms_pfEKF(2,:))
% xlabel('\sigma_{x,y,z} [km/s^2]')
% ylabel('Range Rate Post-fit Residual RMS [km/s]')
% sgtitle('Post-fit Residual Error RMS Values vs \sigma{x,y,z} for an EKF')

% Best case
sigma_xyz_DMC = 1e-9;
Q_DMC = diag([sigma_xyz_DMC, sigma_xyz_DMC, sigma_xyz_DMC].^2);

% Initial priori covariance and state error
Pbar0DMC = diag([1,1,1,1e-3,1e-3,1e-3 1e-16 1e-16 1e-16].^2);
xbar0DMC = [1 1 1 1e-3 1e-3 1e-3 1e-12 1e-12 1e-12]';

% LKF
[XhatLKF_DMC,PLKF_DMC,yLKF_DMC,yhatLKF_DMC] = orbitLKF_DMC(t,xbar0DMC,Pbar0DMC,Q_DMC,noisy_measurements,R,XnomDMC,constants,stations,num_iterations);
plot_filter_diagnostics(t,XnomDMC,XhatLKF_DMC(:,:,end),PLKF_DMC(:,:,:,end),yLKF_DMC(:,:,end),yhatLKF_DMC(:,:,end),R,station_id,"LKF with DMC","")

% EKF
[XhatEKF_DMC,PEKF_DMC,yEKF_DMC,yhatEKF_DMC] = orbitEKF_DMC(t,xbar0DMC,Pbar0DMC,Q_DMC,noisy_measurements,R,XnomDMC,constants,stations,LKFinit);
plot_filter_diagnostics(t,XnomDMC,XhatEKF_DMC(:,:,end),PEKF_DMC(:,:,:,end),yEKF_DMC(:,:,end),yhatEKF_DMC(:,:,end),R,station_id,"EKF with DMC","")

% Acceleration errors
plot_w_diagnostics(t,aJ3,XhatLKF_DMC,PLKF_DMC, "LKF");
plot_w_diagnostics(t,aJ3,XhatEKF_DMC,PEKF_DMC, "EKF");
