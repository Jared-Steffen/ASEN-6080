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
num_iterations = 1;
[XhatLKF1,xhatLKF1,PLKF1,PbarLKF1,yLKF1,yhatLKF1,PhiLKF1] = orbitLKF(t,xbar0,Pbar0,[],"ECI",noisy_measurements,R,Xnom,constants,stations,num_iterations);
plot_filter_diagnostics(t,Xnom,XhatLKF1(:,:,end),PLKF1(:,:,:,end),yLKF1(:,:,end),yhatLKF1(:,:,end),R,station_id,"LKF","")

% SRIF No Process Noise
V = chol(R);
[XhatSRIF1,PSRIF1,ySRIF1,yhatSRIF1] = SRIF(t,xbar0,Pbar0,[],noisy_measurements,V,Xnom,constants,stations);
plot_filter_diagnostics(t,Xnom,XhatSRIF1,PSRIF1,ySRIF1,yhatSRIF1,R,station_id,"SRIF","")

% Compare RMS
Xhat_list{1} = XhatLKF1(:,:,end);
Xhat_list{2} = XhatSRIF1;
plot_compare_RMS(t,Xnom,Xhat_list,{'LKF','SRIF'})

%% Process noise