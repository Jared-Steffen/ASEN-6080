clc; clear; close all

%% Problem 1

% Station lat/lon/alt coordinates
s1_lla = [-35.398333, 148.981944, 0]; % deg, deg, m
s2_lla = [40.427222, 355.739444, 0]; % deg, deg, m
s3_lla = [35.247164, 243.205, 0]; % deg, deg, m
sall_lla = [s1_lla', s2_lla', s3_lla'];
station_ids = 1:length(sall_lla);

% IC
theta0 = deg2rad(122); % deg -> rad

% Earth rotation rate
wE = (2*pi)/24 * 1/3600; % rad/s

% Elevation mask
el_mask = deg2rad(10); % deg -> rad

% genMeasurements struct
measurement_params = struct();
measurement_params.sall_lla    = [s1_lla', s2_lla', s3_lla'];
measurement_params.theta0 = deg2rad(122);          % rad
measurement_params.wE     = (2*pi)/24/3600;        % rad/s
measurement_params.el_mask = deg2rad(10);          % rad

% Grav parameters
mu = 398600.4415;  % km^3/s^2
J2 = 1.0826269e-3;
Re = 6378; % km

% Load truth simulation data - J2
dataJ2 = load("simulation_data.mat");
XnomJ2 = dataJ2.state_unpert;
XnomJ2 = XnomJ2(:,1:6);
XtruthJ2 = dataJ2.state_pert;
XtruthJ2 = XtruthJ2(:,1:6);
PhiJ2 = dataJ2.stm_p2;
PhiJ2 = PhiJ2(1:6,1:6,:);
truth_measurementsJ2 = dataJ2.measurements_real;
nom_measurementsJ2 = dataJ2.measurements_nom;
tJ2 = dataJ2.t;
dx0J2 = dataJ2.pert;
gs_stateJ2 = dataJ2.gs_state;

% Set seed for random number generator
rng(42)

% Measurement covariance matrix
sigma_r = 1e-3; % [km]
sigma_rdot = 1e-6; % [km/s] 
R = diag([sigma_r,sigma_rdot].^2);

% Add noise to measurements
noisy_truth_measurements = cell(length(tJ2),1);
for i = 1:length(tJ2)
    if isempty(truth_measurementsJ2{i})
        continue
    end
    M = truth_measurementsJ2{i};
    n = size(M,1);
    M(2) = M(2) + sigma_r*randn(n,1);
    M(3) = M(3) + sigma_rdot*randn(n,1);
    noisy_truth_measurements{i} = M;
end

% Initial priori covariance and state error
Pbar0 = diag([1,1,1,1e-3,1e-3,1e-3].^2);
xbar0 = dx0J2(1:6);

% Only half the measurements
% half_len = length(noisy_truth_measurements)/2;
% t = t(1:half_len);
% Xnom = Xnom(1:half_len,:);
% Xtruth = Xtruth(1:half_len,:);
% gs_state = gs_state(1:half_len,:);
% nom_measurements = nom_measurements(1:half_len);
% noisy_truth_measurements = noisy_truth_measurements(1:half_len);

% LKF
[dxhat,XhatLKF,PLKF,yLKF,yhatLKF] = orbitLKF(tJ2,xbar0,Pbar0,PhiJ2,noisy_truth_measurements,nom_measurementsJ2,R,XnomJ2,gs_stateJ2);
plot_filter_diagnostics(tJ2,XtruthJ2,XhatLKF,PLKF,yhatLKF,'LKF');

% EKF 
[XhatEKF,PEKF,yEKF,yhatEKF] = orbitEKF(tJ2,xbar0,Pbar0,PhiJ2,noisy_truth_measurements,nom_measurementsJ2,R,XnomJ2,gs_stateJ2,mu,Re,J2,measurement_params);
plot_filter_diagnostics(tJ2,XtruthJ2,XhatEKF,PEKF,yhatEKF,'EKF');

% Batch
tol = 1e-6;
[XhatBLLS,PBLLS,yhatBLLS,batch_cnt] = orbitBatch(tJ2,xbar0,Pbar0,noisy_truth_measurements,R,XnomJ2,gs_stateJ2,mu,Re,J2,measurement_params,tol);
test = XhatBLLS(:,:,end);
% plot_filter_diagnostics(t,Xtruth,XhatBLLS(:,:,1),PBLLS(:,:,:,1),yhatBLLS(:,:,1),'First Batch LLS');
% plot_filter_diagnostics(t,Xtruth,XhatBLLS(:,:,2),PBLLS(:,:,:,2),yhatBLLS(:,:,2),'Second Batch LLS');
% plot_filter_diagnostics(t,Xtruth,XhatBLLS(:,:,3),PBLLS(:,:,:,3),yhatBLLS(:,:,3),'Third Batch LLS');
% plot_filter_diagnostics(t,Xtruth,XhatBLLS(:,:,4),PBLLS(:,:,:,4),yhatBLLS(:,:,4),'Fourth Batch LLS');
plot_filter_diagnostics(tJ2,XtruthJ2,XhatBLLS(:,:,end),PBLLS(:,:,:,end),yhatBLLS(:,:,end),'Final Batch LLS');