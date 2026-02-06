clc; clear; close all

%% Simulation Setup

% Define Earth Grav Model
constants.mu = 3.986004415e5; % km^3/s^2
constants.J2 = 1.0826269e-3;
constants.RE = 6378.1363; % km
constants.wE = 7.2921158553e-5; % rad/s

% Define Atmospheric Drag Model
constants.rho0 = 3.614e-4; % kg/km^3
constants.r0 = 700.0 + constants.RE; % km
constants.H = 88.667; % km
constants.A = 3e-6; % km^2 
constants.m = 970;  % kg
constants.CD = 2.0;

% Define Station Locations (ECEF)
stations.Rs = [-5127.5100, -3794.1600, 0.0; ...
                3860.9100, 3238.4900, 3898.0940; ...
                549.5050, -1380.8720, 6182.1970]; % km
stations.station_ids = [101 337 394]';
stations.el_mask = [];

% Filter Initialization
R = diag([1e-5, 1e-6].^2); % km, km/s
Xr0 = [757.7000 5222.6070 4851.5000]'; % km
Xv0 = [2.21321 4.67834 -5.37130]';  % km/s
X0 = [Xr0; Xv0; constants.mu; constants.J2; constants.CD; stations.Rs(1,:)'; stations.Rs(2,:)'; stations.Rs(3,:)'];

fixed_station = 1; % 1 for station 101, 2 for station 336, 3 for station 394, anything else for none
if fixed_station == 1
    P0 = diag([ ...
        1, 1, 1, ... % position (km^2)
        1, 1, 1, ... % velocity ((km/s)^2)
        1e11, ... % mu (km^3/s^2)^2 
        1e6, ...  % J2
        1e6, ... % CD
        1e-16, 1e-16, 1e-16, ... % station 101 (km^2)
        1, 1, 1, ...  % station 337 (km^2)
        1, 1, 1]); % station 394 (km^2)
elseif fixed_station == 2
    P0 = diag([ ...
        1, 1, 1, ... % position (km^2)
        1, 1, 1, ... % velocity ((km/s)^2)
        1e11, ... % mu (km^3/s^2)^2 
        1e6, ...  % J2
        1e6, ... % CD
        1, 1, 1, ... % station 101 (km^2)
        1e-16, 1e-16, 1e-16, ...  % station 337 (km^2)
        1, 1, 1]); % station 394 (km^2)
elseif fixed_station == 3
    P0 = diag([ ...
        1, 1, 1, ... % position (km^2)
        1, 1, 1, ... % velocity ((km/s)^2)
        1e11, ... % mu (km^3/s^2)^2 
        1e6, ...  % J2
        1e6, ... % CD
        1, 1, 1, ... % station 101 (km^2)
        1, 1, 1, ...  % station 337 (km^2)
        1e-16, 1e-16, 1e-16]); % station 394 (km^2)
else
    P0 = diag([ ...
        1, 1, 1, ... % position (km^2)
        1, 1, 1, ... % velocity ((km/s)^2)
        1e-16, ... % mu (km^3/s^2)^2 
        1e6, ...  % J2
        1e6, ... % CD
        1, 1, 1, ... % station 101 (km^2)
        1, 1, 1, ...  % station 337 (km^2)
        1, 1, 1]); % station 394 (km^2)
end
dx0 = zeros(length(X0),1);

% Load Observation Data
Y = load("project.txt");
Y(:,3:4) = Y(:,3:4)./1000;
t = Y(:,1);
station_id_times = Y(:,2);

% Simulate Initial Orbit
options = odeset('RelTol',1e-11,'AbsTol',1e-11);
[t,Xnom] = ode45(@(t,x) orbitEOM_J2_Drag(t,x,constants),t,X0,options);

%% Batch Filter
tol = 1e-3;
data_types = {"range", "range rate", ""};
data_type = data_types{3}; % 1 for range, 2 for range rate, 3 for both
[XhatBLLS,PBLLS,yBLLS,yhatBLLS,batch_cnt] = orbitBatch(t,dx0,P0,Y,R,Xnom,constants,stations,tol,data_type);

% Simulate Iteration Nominal Orbits
constantsBLLS1.mu = XhatBLLS(1,7,1);
constantsBLLS1.J2 = XhatBLLS(1,8,1);
constantsBLLS1.CD = XhatBLLS(1,9,1);
constantsBLLS1.rho0 = 3.614e-4; % kg/km^3
constantsBLLS1.r0 = 700.0 + constants.RE; % km
constantsBLLS1.H = 88.667; % km
constantsBLLS1.A = 3e-6; % km^2 
constantsBLLS1.m = 970;  % kg
constantsBLLS1.RE = 6378.1363; % km
constantsBLLS1.wE = 7.2921158553e-5; % rad/s
[~,Xnom1] = ode45(@(t,x) orbitEOM_J2_Drag(t,x,constantsBLLS1),t,XhatBLLS(1,:,1),options);
constantsBLLS2.mu = XhatBLLS(1,7,2);
constantsBLLS2.J2 = XhatBLLS(1,8,2);
constantsBLLS2.CD = XhatBLLS(1,9,2);
constantsBLLS2.rho0 = 3.614e-4; % kg/km^3
constantsBLLS2.r0 = 700.0 + constants.RE; % km
constantsBLLS2.H = 88.667; % km
constantsBLLS2.A = 3e-6; % km^2 
constantsBLLS2.m = 970;  % kg
constantsBLLS2.RE = 6378.1363; % km
constantsBLLS2.wE = 7.2921158553e-5; % rad/s
[~,Xnom2] = ode45(@(t,x) orbitEOM_J2_Drag(t,x,constantsBLLS2),t,XhatBLLS(1,:,2),options);
constantsBLLS3.mu = XhatBLLS(1,7,3);
constantsBLLS3.J2 = XhatBLLS(1,8,3);
constantsBLLS3.CD = XhatBLLS(1,9,3);
constantsBLLS3.rho0 = 3.614e-4; % kg/km^3
constantsBLLS3.r0 = 700.0 + constants.RE; % km
constantsBLLS3.H = 88.667; % km
constantsBLLS3.A = 3e-6; % km^2 
constantsBLLS3.m = 970;  % kg
constantsBLLS3.RE = 6378.1363; % km
constantsBLLS3.wE = 7.2921158553e-5; % rad/s
[~,Xnom3] = ode45(@(t,x) orbitEOM_J2_Drag(t,x,constantsBLLS3),t,XhatBLLS(1,:,3),options);

% Plot Diagnostics
plot_filter_diagnostics(t,Xnom1,XhatBLLS(:,:,1),PBLLS(:,:,:,1),yBLLS(:,:,1),yhatBLLS(:,:,1),station_id_times,'First BLLS',data_type)
plot_filter_diagnostics(t,Xnom2,XhatBLLS(:,:,2),PBLLS(:,:,:,2),yBLLS(:,:,2),yhatBLLS(:,:,2),station_id_times,'Second BLLS',data_type)
plot_filter_diagnostics(t,Xnom3,XhatBLLS(:,:,3),PBLLS(:,:,:,3),yBLLS(:,:,3),yhatBLLS(:,:,3),station_id_times,'Third BLLS',data_type)
plot_covariance_ellipsoids(PBLLS(:,:,:,3),'Third BLLS')
plot_covariance_trace(t,PBLLS(:,:,:,3),'Third BLLS');

%% LKF
num_iterations = 1;
[XhatLKF,PLKF,yLKF,yhatLKF] = orbitLKF(t,dx0,P0,Y,R,Xnom,constants,stations,num_iterations);

% Plot Diagnostics
plot_filter_diagnostics(t,Xnom,XhatLKF(:,:,1),PLKF(:,:,:,1),yLKF(:,:,1),yhatLKF(:,:,1),station_id_times,'First LKF',data_type)
plot_covariance_ellipsoids(PLKF(:,:,:,1),'LKF')
plot_covariance_trace(t,PLKF(:,:,:,1),'LKF');
