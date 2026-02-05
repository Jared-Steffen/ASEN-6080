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
P0 = diag([ ...
    1, 1, 1, ... % position (km^2)
    1, 1, 1, ... % velocity ((km/s)^2)
    1e-16, ... % mu (km^3/s^2)^2 
    1e6, ...  % J2
    1e6, ... % CD
    1e-10, 1e-10, 1e-10, ... % station 101 (km^2)
    1, 1, 1, ...  % station 337 (km^2)
    1, 1, 1]); % station 394 (km^2)
dx0 = zeros(length(X0),1);

% Load Observation Data
Y = load("project.txt");
Y(:,3:4) = Y(:,3:4)./1000;
t = Y(:,1);

% Simulate Initial Orbit
options = odeset('RelTol',1e-11,'AbsTol',1e-11);
[t,Xnom1] = ode45(@(t,x) orbitEOM_J2_Drag(t,x,constants),t,X0,options);

%% Batch Filter
tol = 1e-3;
[Xhat,P,y,yhat,batch_cnt] = orbitBatch(t,dx0,P0,Y,R,Xnom1,constants,stations,tol);

% Simulate Iteration Nominal Orbits
[t,Xnom2] = ode45(@(t,x) orbitEOM_J2_Drag(t,x,constants),t,Xhat(1,:,2),options);
[t,Xnom3] = ode45(@(t,x) orbitEOM_J2_Drag(t,x,constants),t,Xhat(1,:,3),options);

plot_filter_diagnostics(t,Xnom1,Xhat(:,:,1),P(:,:,:,1),yhat(:,:,1),'First BLLS')
plot_filter_diagnostics(t,Xnom2,Xhat(:,:,2),P(:,:,:,2),yhat(:,:,2),'Second BLLS')
plot_filter_diagnostics(t,Xnom3,Xhat(:,:,3),P(:,:,:,3),yhat(:,:,3),'Third BLLS')