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

% Load truth simulation data
data = load("simulation_data.mat");
Xnom = data.state_unpert;
Xnom = Xnom(:,1:6);
Xtruth = data.state_pert;
Xtruth = Xtruth(:,1:6);
Phi = data.stm_p2;
Phi = Phi(1:6,1:6,:);
truth_measurements = data.measurements_real;
nom_measurements = data.measurements_nom;
t = data.t;
x0 = data.var;
dx0 = data.pert;
gs_state = data.gs_state;

% Set seed for random number generator
rng(42)

% Measurement covariance matrix
sigma_r = 1e-3; % [km]
sigma_rdot = 1e-6; % [km/s] 
R = diag([sigma_r,sigma_rdot].^2);

% Add noise to measurements
noisy_truth_measurements = cell(length(t),1);
for i = 1:length(t)
    if isempty(truth_measurements{i})
        continue
    end
    M = truth_measurements{i};
    n = size(M,1);
    M(2) = M(2) + sigma_r*randn(n,1);
    M(3) = M(3) + sigma_rdot*randn(n,1);
    noisy_truth_measurements{i} = M;
end

% Initial priori covariance and state error
Pbar0 = diag([1,1,1,1e-3,1e-3,1e-3].^2);
xbar0 = dx0(1:6);

% LKF
% [dxhat,XhatLKF,PLKF,yLKF,yhatLKF] = orbitLKF(t,xbar0,Pbar0,Phi,noisy_truth_measurements,nom_measurements,R,Xnom,gs_state);
% plot_filter_diagnostics(t,Xtruth,XhatLKF,PLKF,yhatLKF,'LKF');

% EKF 
[XhatEKF,PEKF,yEKF,yhatEKF] = orbitEKF(t,xbar0,Pbar0,Phi,noisy_truth_measurements,nom_measurements,R,Xnom,gs_state,mu,Re,J2,measurement_params);
plot_filter_diagnostics(t,Xtruth,XhatEKF,PEKF,yhatEKF,'EKF');

% Batch
tol = 1e-6;
[XhatBLLS,PBLLS,yhatBLLS,batch_cnt] = orbitBatch(t,xbar0,Pbar0,noisy_truth_measurements,R,Xnom,gs_state,mu,Re,J2,measurement_params,tol);
plot_filter_diagnostics(t,Xtruth,XhatBLLS(:,:,4),PBLLS,yhatBLLS(:,:,4),'Batch LLS');