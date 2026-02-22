clc; clear; close all

%% Load data

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

%% Without Process Noise

% Batch
% tol = 0.2;
num_iterations = 1;
% [XhatBLLS,PBLLS,yBLLS,yhatBLLS,batch_cnt] = orbitBatch(t,xbar0,Pbar0,noisy_measurements,R,Xnom,constants,stations,tol,num_iterations,"");
% % plot_filter_diagnostics(t,Xnom,XhatBLLS(:,:,end),PBLLS(:,:,:,end),yBLLS(:,:,end),yhatBLLS(:,:,end),R,station_id,'Batch',"")
% 
% % LKF
% Q = zeros(3);
% [XhatLKF1,xhatLKF1,PLKF1,PbarLKF1,yLKF1,yhatLKF1,PhiLKF1] = orbitLKF(t,xbar0,Pbar0,Q,"ECI",noisy_measurements,R,Xnom,constants,stations,num_iterations);
% % plot_filter_diagnostics(t,Xnom,XhatLKF1(:,:,end),PLKF1(:,:,:,end),yLKF1(:,:,end),yhatLKF1(:,:,end),R,station_id,"LKF","")
% [Xhatl,Pl] = orbitLKFSmoother(t,xhatLKF1,Xnom,PLKF1,PbarLKF1,PhiLKF1);
% % plot_smoother_diagnostics(t,Xnom,Xhatl,Pl)
% 
% % Compare
% Xhat_list{1} = XhatBLLS;
% Xhat_list{2} = XhatLKF1;
% plot_compare_RMS(t,Xnom,Xhat_list,{'Batch (1 Iteration)','Smoothed LKF (No Process Noise)'})

%% With Process Noise

% Best Q
sigma_xyz = 1e-5;
Q = diag([sigma_xyz, sigma_xyz, sigma_xyz].^2);

% LKF
[XhatLKF,xhatLKF,PLKF,PbarLKF,yLKF,yhatLKF,PhiLKF] = orbitLKF(t,xbar0,Pbar0,Q,"ECI",noisy_measurements,R,Xnom,constants,stations,num_iterations);
plot_filter_diagnostics(t,Xnom,XhatLKF(:,:,end),PLKF(:,:,:,end),yLKF(:,:,end),yhatLKF(:,:,end),R,station_id,"LKF","")
[Xhatl,Pl] = orbitLKFSmoother(t,xhatLKF,Xnom,PLKF,PbarLKF,PhiLKF);
plot_smoother_diagnostics(t,Xnom,Xhatl,Pl)
Xhat_list2{1} = Xhatl;
plot_compare_RMS(t,Xnom,Xhat_list2,{'Smoothed LKF'})

% Compare filter to smoother
Xhat_list3{1} = XhatLKF;
Xhat_list3{2} = Xhatl;
plot_compare_RMS(t,Xnom,Xhat_list3,{'LKF','Smoothed LKF'})
plot_covariance_trace(t,PLKF,'LKF')
plot_covariance_trace(t,Pl,'Smoothed LKF')
