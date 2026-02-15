clc; clear; close all

%% Problem 1

% Read in data
file_p1 = fileread('prob1c_solution.json');
data_p1 = jsondecode(file_p1);
% disp(data_p1)

% Constants
mu = data_p1.inputs.state.mu;
J2 = data_p1.inputs.state.J2;
J3 = data_p1.inputs.state.J3;
r_vec = data_p1.inputs.state.r;
Re = 6378; % km

% Outputs
A_ans = data_p1.outputs.A_matrix.values;

% STM fxn call
stm_p1 = STM_J2_J3(mu,J2,J3,Re,r_vec);

%% Problem 2

% Read in data
file_p2 = fileread('prob2b_solution.json');
data_p2 = jsondecode(file_p2);
% disp(data_p2)
truth_2a = load("HW1_truth.txt");

test_X0 = data_p2.inputs.X0.values;
test_Phi0 = data_p2.inputs.Phi0.values;


% Define orbit elements
a = 10000; % km
e = 0.001;
i = deg2rad(40); % deg -> rad
Omega = deg2rad(80); % deg -> rad
w = deg2rad(40); % deg -> rad
nu0 = 0; % rad

% Grav parameters
mu = 398600.4415;  % km^3/s^2
J2 = 1.0826269e-3;
J3 = -2.5324e-6;

% Get initial r and v
[r0,v0] = oe2rv(mu,a,e,Omega,i,w,nu0);
var = [r0;v0;J2];
pert = [1 1 1 1e-3 1e-3 1e-3 0]';
% pert = [0.5 -0.7 0.2 0.0004 -0.0006 0.0002 0]';
var_pert = var + pert;

% Orbital period
n = sqrt(mu/a^3);
T = (2*pi)/n;

% Time vector
tspan = [0 15*T];

% ode45 calls for truth data
options = odeset('RelTol',1e-11,'AbsTol',1e-11);
[t,state_unpert] = ode45(@(t,x) orbitEOM_J2(t,x,mu,Re,J2),tspan,var,options);
[~,state_pert] = ode45(@(t,x) orbitEOM_J2(t,x,mu,Re,J2),t,var_pert,options);

% Reshape ICs for STM
% var_stm = [test_X0;reshape(test_Phi0,[],1)];
var_stm = [var;reshape(eye(7),[],1)];


% ode45 calls for STM
[~,stm] = ode45(@(t,x) odeSTM_J2(t,x,mu,Re),t,var_stm,options);

% Extract STMs and reshape to be 7x7 and then propogate pertubations
for i = 1:length(t)
    stm_p2(:,:,i) = reshape(stm(i,8:end),7,7);
    pertt(i,:) = stm_p2(:,:,i) * pert;
end


% Plots
% plot_rv_state(t,test_states,{"Given Data","Simulated Data"},false)
% plot_rv_state(t,deltaxs,{"NL Propagation","STM Propagation"},true)

r_labels = {'\deltax Position Difference [km]','\deltay Position Difference [km]','\deltaz Position Difference [km]'};
v_labels = {'\deltav_x Velocity Difference [km/s]','\deltav_y Velocity Difference [km/s]','\deltav_z Velocity Difference [km/s]'};

% pert_diff = (state_pert - state_unpert) - pertt;
% r = pert_diff(:,1:3);
% v = pert_diff(:,4:6);

% figure();
% subplot(311)
% hold on
% plot(t/3600,r(:,1))
% xlabel('Time [hr]')
% ylabel(r_labels{1})
% subplot(312)
% hold on
% plot(t/3600,r(:,2))
% xlabel('Time [hr]')
% ylabel(r_labels{2})
% subplot(313)
% hold on
% plot(t/3600,r(:,3))
% xlabel('Time [hr]')
% ylabel(r_labels{3})
% sgtitle('\deltar_{NL} - \deltar_{STM} Position Difference')
% 
% figure();
% subplot(311)
% hold on
% plot(t/3600,v(:,1))
% xlabel('Time [hr]')
% ylabel(v_labels{1})
% subplot(312)
% hold on
% plot(t/3600,v(:,2))
% xlabel('Time [hr]')
% ylabel(v_labels{2})
% subplot(313)
% hold on
% plot(t/3600,v(:,3))
% xlabel('Time [hr]')
% ylabel(v_labels{3})
% sgtitle('\deltav_{NL} - \deltav_{STM} Velocity Difference')

%% Problem 3

% Read in data
file_p3b = fileread('prob3b_solution.json');
data_p3b = jsondecode(file_p3b);
disp(data_p3b)

test_sc_state = [data_p3b.inputs.spacecraft_state.r;data_p3b.inputs.spacecraft_state.v];
test_gs_state = [data_p3b.inputs.station_state.Rs;data_p3b.inputs.station_state.Vs];

% Verify fxn
test_sc_Htilde = sc_range_ranger_Htilde(test_sc_state,test_gs_state);

% Read in data
file_p3d = fileread('prob3d_solution.json');
data_p3d = jsondecode(file_p3d);
disp(data_p3d)

% Earth rotation rate
wE = (2*pi)/24 * 1/3600; % rad/s

% Verify fxn
test_gs_Htilde = gs_range_ranger_Htilde(test_sc_state,test_gs_state,wE);


%% Problem 4

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

% Generate measurements
[measurements_real,gs_state] = genMeasurements(t,sall_lla,theta0,wE,el_mask,state_pert);
[measurements_nom,~] = genMeasurements(t,sall_lla,theta0,wE,el_mask,state_unpert);
measurements_shift = measurements_real; % copy for adding Doppler calculations
measurements_noisy = measurements_real; % copy for adding noise

% Save measurements for HW 2
save('simulation_dataJ2_test.mat','t','state_pert','state_unpert','measurements_real','measurements_nom','var','pert',"gs_state")

% Plot a and b
plot_measurements(t,measurements_real,station_ids)

% Doppler and RU conversion constants
fTref = 8.44e9; % Hz
c = 3e8; % m/s

% Extract range and range rate measurements
for i = 1:length(measurements_real)
    M = measurements_shift{i};
    if isempty(M)
        continue
    end
    M(:,2) = 221/749*(M(:,2).*1000)/c*fTref;
    M(:,3) = (-2*(M(:,3).*1000))/c*fTref;
    measurements_shift{i} = M;
end

% Plot c

% Generate legend
station_legend = cell(length(station_ids),1);
for i = 1:length(station_ids)
    station_legend{i} = "Station ID: "+ station_ids(i);
end

% Initialize
Ns = 3;
T  = numel(t);
station_id_all = [];
time_all = [];
RU_all = [];
fshift_all = [];
el_all = [];

for i = 1:T
    M = measurements_shift{i};
    if isempty(M) 
        continue; 
    end
    n = size(M,1);
    station_id_all = [station_id_all; M(:,1)];
    time_all = [time_all; repmat(t(i), n, 1)];
    RU_all = [RU_all; M(:,2)];
    fshift_all = [fshift_all; M(:,3)];
end

% Range and range rate figure
figure();
subplot(211) 
hold on
for s = 1:Ns
    idx = (station_id_all == s);
    if any(idx)
        plot(time_all(idx)/3600, RU_all(idx), 'o', 'MarkerSize', 4);
    end
end
ylabel('Range Units')
xlabel('Time [hr]')
title('Range')

subplot(212)
hold on
for s = 1:Ns
    idx = (station_id_all == s);
    if any(idx)
        plot(time_all(idx)/3600, fshift_all(idx), 'o', 'MarkerSize', 4);
    end
end
ylabel('Doppler Shift [Hz]')
xlabel('Time [hr]')
title('Range Rate')
lgd = legend(station_legend);
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

% Noise
sigma_rhodot = 0.5e-6; % km/s

% Extract range rate and add noise to measurements
noisy_residuals = cell(T,1);
for i = 1:length(measurements_noisy)
    M = measurements_noisy{i};
    Mn = measurements_real{i};
    if isempty(M)
        continue
    end
    n = size(M,1);
    noise = sigma_rhodot*randn(n,1);
    M(:,3) = M(:,3)+noise;
    measurements_noisy{i} = M;
    noisy_residuals{i} = [M(:,1),M(:,3) - Mn(:,3)];
end

% Plot noisy measurements
plot_measurements(t,measurements_noisy,station_ids)

% Initialize
Ns = 3;
T  = numel(t);
station_id_all = [];
time_all = [];
rho_dot_resid_all = [];

for i = 1:T
    M = noisy_residuals{i};
    if isempty(M) 
        continue; 
    end
    n = size(M,1);
    station_id_all = [station_id_all; M(:,1)];
    time_all = [time_all; repmat(t(i), n, 1)];
    rho_dot_resid_all = [rho_dot_resid_all; M(:,2)];
end

figure();
hold on
for s = 1:Ns
    idx = (station_id_all == s);
    if any(idx)
        plot(time_all(idx)/3600, rho_dot_resid_all(idx), 'o', 'MarkerSize', 4);
    end
end
ylabel('Range Rate Error [km/s]')
xlabel('Time [hr]')
title('Range Rate Noise Difference')
lgd = legend(station_legend);
lgd.Units = 'normalized';
lgd.Position = [0.75 0.935 0.25 0.05];

% For HW 2
[tJ3,state_unpertJ3] = ode45(@(t,x) orbitEOM_J2_J3(t,x,mu,Re,J2,J3),tspan,var,options);
[~,state_pertJ3] = ode45(@(t,x) orbitEOM_J2_J3(t,x,mu,Re,J2,J3),tJ3,var_pert,options);
[measurements_realJ3,gs_stateJ3] = genMeasurements(tJ3,sall_lla,theta0,wE,el_mask,state_pertJ3);
[measurements_nomJ3,~] = genMeasurements(tJ3,sall_lla,theta0,wE,el_mask,state_unpertJ3);

% Save measurements for HW 2
save('simulation_dataJ3.mat','tJ3','state_pertJ3','state_unpertJ3','measurements_realJ3','measurements_nomJ3','var','pert',"gs_stateJ3")

% Plot differences
meas_diff = cell(size(measurements_realJ3));
for k = 1:numel(meas_diff)
    mr = measurements_real{k};
    mrJ3 = measurements_realJ3{k};
    if isempty(mr)
        continue
    end
    meas_diff{k} = [mr(1), mr(2:end) - mrJ3(2:end)];
end
plot_measurements(tJ3,meas_diff,station_ids)

diffJ3 = {state_pert - state_pertJ3};
plot_rv_state(tJ3,diffJ3,[],true)