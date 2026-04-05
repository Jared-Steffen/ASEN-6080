clc; clear; close all

%% Simulation Setup
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% Time
t = 0:60:3600*24; % 24 hours

% Define Earth Grav Model
constants.mu = 3.986004415e5; % km^3/s^2
constants.J2 = 1.0826269e-3;
constants.RE = 6378.1363; % km
constants.wE = 7.2921158553e-5; % rad/s

% Define Atmospheric Drag Model
constants.rho0 = 3.614e-4; % kg/km^3
constants.r0 = 500.0 + constants.RE; % km
constants.H = 88.667; % km
constants.A = 3e-6; % km^2 
constants.m = 970;  % kg
constants.CD = 2.0;

% Initial Conditions
v0 = sqrt(constants.mu / constants.r0); % circular velocity
X0 = [constants.r0; 0; 0; 0; v0*cosd(45); v0*sind(45)];

% Covariance
sigma_r = 1; % km
sigma_v = 1e-3; % km/s

P0 = diag([sigma_r^2 sigma_r^2 sigma_r^2 ...
           sigma_v^2 sigma_v^2 sigma_v^2 ]);
%% Monte Carlo Simulation
n = 6;
N = 1000;
for i = 1:N
    X0_pert = X0 + chol(P0,'lower') * randn(6,1);
    [~,X(:,:,i)] = ode45(@(t,x) odeSTM_J2_Drag(t,x,constants),t,[X0_pert;reshape(eye(n),[],1)],options); 
end

    
%% Nominal Simulation 
[~,Xnom] = ode45(@(t,x) odeSTM_J2_Drag(t,x,constants),t,[X0;reshape(eye(n),[],1)],options); 

%% Corner Plots
labels = {'x','y','z','vx','vy','vz'};
t_find = [0,6,12,18,24]*3600;
t_idxs = zeros(size(t_find));

Xmc = X(:,1:6,:);

for i = 1:length(t_find)
    t_idx = find(t == t_find(i));
    cornerPlot(Xmc,Xnom(:,1:6),t,t_idx,labels)
end

%% Linear Covariance Propagation

for i = 1:length(t)
    Phi = reshape(Xnom(i,7:end),6,6);
    Pk_LKF(:,:,i) = Phi*P0*Phi';
    Xk_LKF(:,i) = Phi*Xnom(i,1:6)';
end

for i = 1:length(t_find)
    t_idx = find(t == t_find(i));
    cornerPlotFilter(Xmc,Xnom(:,1:6),Pk_LKF,t,t_idx,labels)
    twoSigmaStats(Xmc,Xnom(:,1:6),Pk_LKF,t,t_idx,labels)
end


%% UKF Covariance Propagation

% Form constants and weights
alpha = 0.1;
kappa = -3;
beta = 2;
lambda = alpha^2*(n+kappa)-n;
gamma = sqrt(n+lambda);
Wm0 = lambda/(n+lambda);
Wc0 = lambda/(n+lambda)+(1-alpha^2+beta);
Wmc = 1/(2*(n+lambda));
Wm = [Wm0;Wmc*ones(2*n,1)];
Wc = [Wc0;Wmc*ones(2*n,1)];

S = chol(P0,'lower');
Chi0 = reshape([X0 X0+gamma.*S X0-gamma.*S],[],1);  
[~,Chii] = ode45(@(t,x) orbitEOM_J2_Drag_UKF(t,x,constants),t,Chi0,options);

for i = 1:length(t)
    % Get new mean
    Chi = reshape(Chii(i,:),n,2*n+1);
    sum_m = zeros(size(X0(:)));
    for j = 1:size(Chi,2)
        sum_m = sum_m + Wm(j).*Chi(:,j); 
    end
    Xk_UKF(:,i) = sum_m;
    
    % Get new covariance
    sum_c = zeros(size(P0));
    for j = 1:size(Chi,2)
        sum_c = sum_c + Wc(j).*(Chi(:,j)-sum_m)*(Chi(:,j)-sum_m)';
    end
    Pk_UKF(:,:,i) = sum_c;
end

for i = 1:length(t_find)
    t_idx = find(t == t_find(i));
    cornerPlotFilter(Xmc,Xnom(:,1:6),Pk_UKF,t,t_idx,labels)
    twoSigmaStats(Xmc,Xnom(:,1:6),Pk_UKF,t,t_idx,labels)
end