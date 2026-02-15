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
n = length(X0);

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
        1e11, ... % mu (km^3/s^2)^2 
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
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,Xnom] = ode45(@(t,x) odeSTM_J2_Drag(t,x,constants),t,[X0;reshape(eye(n),[],1)],options);

%% Batch Filter
tol = 1e-3;
data_types = {"range", "range rate", ""}; % 1 for range, 2 for range rate, 3 for both
% [XhatBLLS,PBLLS,yBLLS,yhatBLLS,batch_cnt] = orbitBatch(t,dx0,P0,Y,R,Xnom,constants,stations,tol,data_types{3});
% 
% % Simulate Iteration Nominal Orbits
% constantsBLLS1.mu = XhatBLLS(1,7,1);
% constantsBLLS1.J2 = XhatBLLS(1,8,1);
% constantsBLLS1.CD = XhatBLLS(1,9,1);
% constantsBLLS1.rho0 = 3.614e-4; % kg/km^3
% constantsBLLS1.r0 = 700.0 + constants.RE; % km
% constantsBLLS1.H = 88.667; % km
% constantsBLLS1.A = 3e-6; % km^2 
% constantsBLLS1.m = 970;  % kg
% constantsBLLS1.RE = 6378.1363; % km
% constantsBLLS1.wE = 7.2921158553e-5; % rad/s
% [~,Xnom1] = ode45(@(t,x) odeSTM_J2_Drag(t,x,constantsBLLS1),t,[XhatBLLS(1,:,1)';reshape(eye(n),[],1)],options);
% constantsBLLS2.mu = XhatBLLS(1,7,2);
% constantsBLLS2.J2 = XhatBLLS(1,8,2);
% constantsBLLS2.CD = XhatBLLS(1,9,2);
% constantsBLLS2.rho0 = 3.614e-4; % kg/km^3
% constantsBLLS2.r0 = 700.0 + constants.RE; % km
% constantsBLLS2.H = 88.667; % km
% constantsBLLS2.A = 3e-6; % km^2 
% constantsBLLS2.m = 970;  % kg
% constantsBLLS2.RE = 6378.1363; % km
% constantsBLLS2.wE = 7.2921158553e-5; % rad/s
% [~,Xnom2] = ode45(@(t,x) odeSTM_J2_Drag(t,x,constantsBLLS2),t,[XhatBLLS(1,:,2)';reshape(eye(n),[],1)],options);
% constantsBLLS3.mu = XhatBLLS(1,7,3);
% constantsBLLS3.J2 = XhatBLLS(1,8,3);
% constantsBLLS3.CD = XhatBLLS(1,9,3);
% constantsBLLS3.rho0 = 3.614e-4; % kg/km^3
% constantsBLLS3.r0 = 700.0 + constants.RE; % km
% constantsBLLS3.H = 88.667; % km
% constantsBLLS3.A = 3e-6; % km^2 
% constantsBLLS3.m = 970;  % kg
% constantsBLLS3.RE = 6378.1363; % km
% constantsBLLS3.wE = 7.2921158553e-5; % rad/s
% [~,Xnom3] = ode45(@(t,x) odeSTM_J2_Drag(t,x,constantsBLLS3),t,[XhatBLLS(1,:,3)';reshape(eye(n),[],1)],options);
% 
% % Plot Diagnostics
% plot_filter_diagnostics(t,Xnom1(:,1:n),XhatBLLS(:,:,1),PBLLS(:,:,:,1),yBLLS(:,:,1),yhatBLLS(:,:,1),R,station_id_times,'First NLLS',data_types{3})
% plot_filter_diagnostics(t,Xnom3(:,1:n),XhatBLLS(:,:,3),PBLLS(:,:,:,3),yBLLS(:,:,3),yhatBLLS(:,:,3),R,station_id_times,'Third NLLS',data_types{3})
% plot_constant_error(t,Xnom3(:,1:n),XhatBLLS(:,:,3),PBLLS(:,:,:,3),'Third NLLS')
% plot_station_error(t,Xnom3(:,1:n),XhatBLLS(:,:,3),PBLLS(:,:,:,3),'Third NLLS',stations.station_ids)
% plot_covariance_trace(t,PBLLS(:,:,:,3),'Third NLLS');
% 
% %% LKF
% num_iterations = 1;
% [XhatLKF,PLKF,yLKF,yhatLKF] = orbitLKF(t,dx0,P0,Y,R,Xnom,constants,stations,num_iterations);
% 
% % Plot Diagnostics
% plot_filter_diagnostics(t,Xnom(:,1:n),XhatLKF(:,:,end),PLKF(:,:,:,end),yLKF(:,:,end),yhatLKF(:,:,end),R,station_id_times,'CKF',data_types{3})
% plot_constant_error(t,Xnom(:,1:n),XhatLKF(:,:,end),PLKF(:,:,:,end),'CKF')
% plot_station_error(t,Xnom(:,1:n),XhatLKF(:,:,end),PLKF(:,:,:,end),'CKF',stations.station_ids)
% plot_covariance_trace(t,PLKF(:,:,:,end),'CKF');
% 
% % Compare Ellipsoids
% plot_covariance_ellipsoids({PBLLS(:,:,:,3),PLKF(:,:,:,1)},{'Third NLLS','CKF'})
% 
% %% Data Strength
% [XhatBLLS2,PBLLS2,yBLLS2,yhatBLLS2,batch_cnt2] = orbitBatch(t,dx0,P0,Y,R,Xnom,constants,stations,tol,data_types{1});
% [XhatBLLS3,PBLLS3,yBLLS3,yhatBLLS3,batch_cnt3] = orbitBatch(t,dx0,P0,Y,R,Xnom,constants,stations,tol,data_types{2});
% 
% % Plot Diagnostics
% constantsBLLS4.mu = XhatBLLS2(1,7,1);
% constantsBLLS4.J2 = XhatBLLS2(1,8,1);
% constantsBLLS4.CD = XhatBLLS2(1,9,1);
% constantsBLLS4.rho0 = 3.614e-4; % kg/km^3
% constantsBLLS4.r0 = 700.0 + constants.RE; % km
% constantsBLLS4.H = 88.667; % km
% constantsBLLS4.A = 3e-6; % km^2 
% constantsBLLS4.m = 970;  % kg
% constantsBLLS4.RE = 6378.1363; % km
% constantsBLLS4.wE = 7.2921158553e-5; % rad/s
% [~,Xnom4] = ode45(@(t,x) odeSTM_J2_Drag(t,x,constantsBLLS4),t,[XhatBLLS2(1,:,1)';reshape(eye(n),[],1)],options);
% constantsBLLS5.mu = XhatBLLS3(1,7,2);
% constantsBLLS5.J2 = XhatBLLS3(1,8,2);
% constantsBLLS5.CD = XhatBLLS3(1,9,2);
% constantsBLLS5.rho0 = 3.614e-4; % kg/km^3
% constantsBLLS5.r0 = 700.0 + constants.RE; % km
% constantsBLLS5.H = 88.667; % km
% constantsBLLS5.A = 3e-6; % km^2 
% constantsBLLS5.m = 970;  % kg
% constantsBLLS5.RE = 6378.1363; % km
% constantsBLLS5.wE = 7.2921158553e-5; % rad/s
% [~,Xnom5] = ode45(@(t,x) odeSTM_J2_Drag(t,x,constantsBLLS5),t,[XhatBLLS3(1,:,2)';reshape(eye(n),[],1)],options);
% 
% plot_filter_diagnostics(t,Xnom4(:,1:n),XhatBLLS2(:,:,3),PBLLS2(:,:,:,3),yBLLS2(:,:,3),yhatBLLS2(:,:,3),R,station_id_times,'Third NLLS - Range Only',data_types{1})
% plot_constant_error(t,Xnom4(:,1:n),XhatBLLS2(:,:,3),PBLLS2(:,:,:,3),'Third NLLS - Range Only')
% plot_station_error(t,Xnom4(:,1:n),XhatBLLS2(:,:,3),PBLLS2(:,:,:,3),'Third NLLS - Range Only',stations.station_ids)
% plot_covariance_trace(t,PBLLS2(:,:,:,3),'Third NLLS - Range Only');
% 
% plot_filter_diagnostics(t,Xnom5(:,1:n),XhatBLLS3(:,:,3),PBLLS3(:,:,:,3),yBLLS3(:,:,3),yhatBLLS3(:,:,3),R,station_id_times,'Third NLLS - Range Rate Only',data_types{2})
% plot_constant_error(t,Xnom5(:,1:n),XhatBLLS3(:,:,3),PBLLS3(:,:,:,3),'Third NLLS - Range Rate Only')
% plot_station_error(t,Xnom5(:,1:n),XhatBLLS3(:,:,3),PBLLS3(:,:,:,3),'Third NLLS - Range Rate Only',stations.station_ids)
% plot_covariance_trace(t,PBLLS3(:,:,:,3),'Third NLLS - Range Rate Only');
% 
% % Compare Ellipsoids
% plot_covariance_ellipsoids({PBLLS(:,:,:,3),PBLLS2(:,:,:,3),PBLLS3(:,:,:,3)},{'Range and Range Rate NLLS','Range NLLS','Range Rate NLLS'})

%% EKF 
LKFinit = 100;
[XhatEKF,PEKF,yEKF,yhatEKF] = orbitEKF(t,dx0,P0,Y,R,Xnom(:,1:n),constants,stations,LKFinit);

% Plot Diagnostics
plot_filter_diagnostics(t,Xnom(:,1:n),XhatEKF,PEKF,yEKF,yhatEKF,R,station_id_times,'EKF',data_types{3})
plot_constant_error(t,Xnom(:,1:n),XhatEKF,PEKF,'EKF')
plot_station_error(t,Xnom(:,1:n),XhatEKF,PEKF,'EKF',stations.station_ids)
plot_covariance_trace(t,PEKF,'EKF');

% Compare Ellipsoids
plot_covariance_ellipsoids({PBLLS(:,:,:,3),PLKF(:,:,:,1),PEKF},{'Third NLLS','CKF','EKF'})