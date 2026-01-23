clc; clear; close all

%% Problem 1

% Load truth simulation data
data = load("simulation_data.mat");
state_nom = data.state_unpert;
state_truth = data.state_pert;
stm_nom = data.stm_p2;
truth_measurements = data.measurements;
t = data.t;
x0 = data.var;
dx0 = data.pert;

% Set seed for random number generator
rng(42)

% Measurement variances
sigma_r = 1e-3; % [km]
sigma_rdot = 1e-6; % [km/s] 

% Add noise to measurements
noisy_truth_measurements = cell(length(t),1);
for i = 1:length(t)
    if isempty(truth_measurements{i})
        continue
    end
    M = truth_measurements{i};
    n = size(M,1);
    M(2) = M(2) + sigma_r*rand(n,1);
    M(3) = M(3) + sigma_rdot*randn(n,1);
    noisy_truth_measurements{i} = M;
end

% Initial priori covariance and state error
Pbar0 = diag([1,1,1,1e-3,1e-3,1e-3]);
xbar0 = diag([0,0,0,0,0,0]);

% LKF


