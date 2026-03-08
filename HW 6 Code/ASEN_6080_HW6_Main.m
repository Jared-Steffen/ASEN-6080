clc; clear; close all

%% Load data

% HW J3 Data
data = load("simulation_dataJ2J3_test.mat");

% Split up data
t = data.t;
Xnom = data.Xnom;
noisy_measurements = data.noisy_measurements;
station_id = noisy_measurements(:,2);
R = data.R;
Pbar0 = data.Pbar0;
xbar0 = data.xbar0;
X0 = Xnom(1,:) + xbar0';
constants = data.constants;
stations = data.stations;


%% No process noise

% UKF
% alpha = 1;
% [XhatUKF,PUKF,yUKF,yhatUKF] = UKF_J2(t,X0,Pbar0,[],alpha,noisy_measurements,R,constants,stations);
% plot_filter_diagnostics(t,Xnom,XhatUKF,PUKF,yUKF,yhatUKF,R,station_id,"UKF","")
% 
% % EKF
% LKFinit = 10;
% [XhatEKF,PEKF,yEKF,yhatEKF] = orbitEKF(t,xbar0,Pbar0,[],"ECI",noisy_measurements,R,Xnom,constants,stations,LKFinit);
% plot_filter_diagnostics(t,Xnom,XhatEKF,PEKF,yEKF,yhatEKF,R,station_id,"EKF","")

%% Process noise

% Process noise
% sigma_xyz = 1e-5;
% Q = diag([sigma_xyz, sigma_xyz, sigma_xyz].^2);
% 
% % UKF 1
% alpha = 1;
% [XhatUKF2,PUKF2,yUKF2,yhatUKF2] = UKF_J2(t,X0,Pbar0,Q,alpha,noisy_measurements,R,constants,stations);
% plot_filter_diagnostics(t,Xnom,XhatUKF2,PUKF2,yUKF2,yhatUKF2,R,station_id,"UKF","")
% 
% % UKF 2
% alpha = 1e-4;
% [XhatUKF3,PUKF3,yUKF3,yhatUKF3] = UKF_J2(t,X0,Pbar0,Q,alpha,noisy_measurements,R,constants,stations);
% plot_filter_diagnostics(t,Xnom,XhatUKF3,PUKF3,yUKF3,yhatUKF3,R,station_id,"UKF","")
% 
% % EKF
% LKFinit = 10;
% [XhatEKF2,PEKF2,yEKF2,yhatEKF2] = orbitEKF(t,xbar0,Pbar0,Q,"ECI",noisy_measurements,R,Xnom,constants,stations,LKFinit);
% plot_filter_diagnostics(t,Xnom,XhatEKF2,PEKF2,yEKF2,yhatEKF2,R,station_id,"EKF","")
% 
% % Comparison
% XhatComp{1} = XhatEKF2;
% XhatComp{2} = XhatUKF2;
% XhatComp{3} = XhatUKF3;
% plot_compare_RMS(t,Xnom,XhatComp,{"EKF","UKF with \alpha = 1", "UKF with \alpha = 1e-4"})

%% UKF with larger initial errors

% Process noise
sigma_xyz = 1e-5;
Q = diag([sigma_xyz, sigma_xyz, sigma_xyz].^2);

% Perturb ICs more
% xbar0 = 10.*xbar0;
% X0 = Xnom(1,:) + xbar0';
% Pbar0 = 100.*Pbar0;

% UKF 1
alpha = 1e-4;
[XhatUKF2,PUKF2,yUKF2,yhatUKF2] = UKF_J2(t,X0,Pbar0,Q,alpha,noisy_measurements,R,constants,stations);
plot_filter_diagnostics(t,Xnom,XhatUKF2,PUKF2,yUKF2,yhatUKF2,R,station_id,"UKF","")

% EKF
LKFinit = 10;
% [XhatEKF2,PEKF2,yEKF2,yhatEKF2] = orbitEKF(t,xbar0,Pbar0,Q,"ECI",noisy_measurements,R,Xnom,constants,stations,LKFinit);
% plot_filter_diagnostics(t,Xnom,XhatEKF2,PEKF2,yEKF2,yhatEKF2,R,station_id,"EKF","")

% Comparison
% XhatComp{1} = XhatEKF2;
% XhatComp{2} = XhatUKF2;
% plot_compare_RMS(t,Xnom,XhatComp,{"EKF","UKF with \alpha = 1e-4"})

%% UKF with J3 dynamics

% UKF
[XhatUKF3,PUKF3,yUKF3,yhatUKF3] = UKF_J3(t,X0,Pbar0,Q,alpha,noisy_measurements,R,constants,stations);
plot_filter_diagnostics(t,Xnom,XhatUKF3,PUKF3,yUKF3,yhatUKF3,R,station_id,"UKF","")

% Comparison
XhatComp{1} = XhatUKF2;
XhatComp{2} = XhatUKF3;
plot_compare_RMS(t,Xnom,XhatComp,{"UKF without J3","UKF with J3"})