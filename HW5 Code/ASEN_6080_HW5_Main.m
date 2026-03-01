clc; clear; close all

%% No process noise

% Load in HW J2 Data
data = load("simulation_dataJ2_test.mat");

% Split up data
t = data.tJ2;
Xnom = data.XnomJ2;
noisy_measurements = data.noisy_measurementsJ2;
station_id = noisy_measurements(:,2);
R = data.R;
Pbar0 = data.Pbar0;
xbar0 = data.xbar0;
constants = data.constants;
stations = data.stations;

% LKF No Process Noise
% num_iterations = 1;
% [XhatLKF1,xhatLKF1,PLKF1,PbarLKF1,yLKF1,yhatLKF1,PhiLKF1] = orbitLKF(t,xbar0,Pbar0,[],"ECI",noisy_measurements,R,Xnom,constants,stations,num_iterations);
% plot_filter_diagnostics(t,Xnom,XhatLKF1(:,:,end),PLKF1(:,:,:,end),yLKF1(:,:,end),yhatLKF1(:,:,end),R,station_id,"LKF","")

% SRIF No Process Noise
% V = chol(R);
% [XhatSRIF1,PSRIF1,ySRIF1,yhatSRIF1] = SRIF(t,xbar0,Pbar0,[],noisy_measurements,V,Xnom,constants,stations);
% plot_filter_diagnostics(t,Xnom,XhatSRIF1,PSRIF1,ySRIF1,yhatSRIF1,R,station_id,"SRIF","")

% Compare RMS
% Xhat_list{1} = XhatLKF1(:,:,end);
% Xhat_list{2} = XhatSRIF1;
% plot_compare_RMS(t,Xnom,Xhat_list,{'LKF','SRIF'})

%% Process noise

% Load in HW J3 Data
data = load("simulation_dataJ2J3_test.mat");

% Split up data
t = data.t;
Xnom = data.Xnom;
noisy_measurements = data.noisy_measurements;
station_id = noisy_measurements(:,2);
R = data.R;
Pbar0 = data.Pbar0;
xbar0 = data.xbar0;
constants = data.constants;
stations = data.stations;

% Best Q
sigma_xyz = 1e-5;
Q = diag([sigma_xyz, sigma_xyz, sigma_xyz].^2);

% LKF w/ Process Noise
num_iterations = 1;
[XhatLKF2,xhatLKF2,PLKF2,PbarLKF2,yLKF2,yhatLKF2,PhiLKF2] = orbitLKF(t,xbar0,Pbar0,Q,"ECI",noisy_measurements,R,Xnom,constants,stations,num_iterations);
plot_filter_diagnostics(t,Xnom,XhatLKF2(:,:,end),PLKF2(:,:,:,end),yLKF2(:,:,end),yhatLKF2(:,:,end),R,station_id,"LKF","")

% SRIF w/ Process Noise
V = chol(R);
[XhatSRIF2,PSRIF2,ySRIF2,yhatSRIF2] = SRIF(t,xbar0,Pbar0,Q,noisy_measurements,V,Xnom,constants,stations);
plot_filter_diagnostics(t,Xnom,XhatSRIF2,PSRIF2,ySRIF2,yhatSRIF2,R,station_id,"SRIF","")

% Compare RMS
% Xhat_list2{1} = XhatLKF2(:,:,end);
% Xhat_list2{2} = XhatSRIF2;
% plot_compare_RMS(t,Xnom,Xhat_list2,{'LKF','SRIF'})