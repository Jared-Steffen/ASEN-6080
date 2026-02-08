clear
clc
close all

% addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")
load("data.txt")
data(:,3:4) = data(:,3:4)./1000; % convert to km

% constants 
ae = 6378136.3e-3; % [km] mean equitorial radius of earth
constants.ae = ae; % [km] mean equitorial radius of earth
area = 3e-6; % km^2
constants.area = area;
rho0 = 3.614e-4; % kg/km^3
constants.rho0 = rho0;
r0 = 700000.0e-3 + ae; % km
constants.r0 = r0;
H0 = 88667.0e-3; % km
constants.H0 = H0;
m = 970; % kg
constants.m = m;
constants.omegaE = 7.2921158553e-5; % [rad/s]
constants.theta0 = 0; % [rad]


%% Simulate True Trajectory

n = 18; % length of state vector
% Initial Conditions
r0_N = [757.700; 5222.607; 4851.500]; % [km]
v0_N = [2.21321; 4.67834; -5.3713]; % [km/s]
mu = 3.986004415e5; % [km^3/s^2] earth's gravitational parameter
J2 = 1.082626925638815e-3; % J2 perturbation
Cd = 2; % approximate coefficient of drag
R0s1_E = [-5127.5100; -3794.1600; 0.0]; % [km]
R0s2_E = [3860.9100; 3238.4900; 3898.0940]; % [km]
R0s3_E = [549.5050; -1380.8720; 6182.1970]; % [km]

X0 = [r0_N; v0_N; mu; J2; Cd; R0s1_E; R0s2_E; R0s3_E];
S0 = [X0; reshape(eye(n),[],1)];

% get time span
t = data(:,1);

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~,Sref] = ode45(@(t,S) dragJ2ODE(t,S,constants), t, S0, options);
% plotEarthOrbit(Sref(:,1:3)', ae, "Simulated Reference Trajectory")
Xref = Sref(:,1:n)';

dx0 = zeros(n,1);
P0 = diag([ones(6,1); 1e11; 1e6; 1e6; 1e-16*ones(3,1); ones(6,1)]); 
noise_sd = [1e-5; 1e-6]; % km km/s
R = [noise_sd(1)^2, 0; 0, noise_sd(2)^2];

max_iterations = 8;

%% CKF

[X_CKF, dx_CKF, P_CKF, y_CKF, alpha_CKF] = newCKF(data, R, X0, dx0, P0, 1, "both", constants);
titles = ["CKF Pre-Fit Residuals vs. Time", "CKF Post-Fit Residuals vs. Time"];
plotResiduals(data(:,1), y_CKF, alpha_CKF, noise_sd, titles)
% RMSresidual_CKF = sqrt(1/length(t)*sum([alpha_CKF(1,:);alpha_CKF(2,:)].^2, 2));
% RMSresidualpre_CKF = sqrt(1/length(t)*sum([y_CKF(1,:);y_CKF(2,:)].^2, 2));
sd_CKF = diag(P_CKF(:,:,end))'.^(1/2);
Xf_CKF = X_CKF(:,end)' + dx_CKF(:,end)';

trueX_CKF = X_CKF+dx_CKF;
deltaX_CKF = Xref -(X_CKF+dx_CKF);
plotDeltaX(t, deltaX_CKF, "CKF $\delta x_{LS}$ vs Time")
% plotDeltaXsc(t, deltaX_CKF, "CKF $\delta x_{LS}$ vs Time")
% 
% [X_CKFiterated, dx_CKFiterated, P_CKFiterated, y_CKFiterated, alpha_CKFiterated] = newCKF(data, R, X0, dx0, P0, 20, "both", constants);
% titles = ["Iterated CKF Pre-Fit Residuals vs. Time", "Iterated CKF Post-Fit Residuals vs. Time"];
% plotResiduals(data(:,1), y_CKFiterated, alpha_CKFiterated, noise_sd, titles)
% RMSresidual_iteratedCKF = sqrt(1/length(t)*sum([alpha_CKFiterated(1,:);alpha_CKFiterated(2,:)].^2, 2));
% RMSresidualpre_iteratedCKF = sqrt(1/length(t)*sum([y_CKFiterated(1,:);y_CKFiterated(2,:)].^2, 2));
% sd_CKFiterated = diag(P_CKFiterated(:,:,end))'.^(1/2);
% Xf_CKFiterated = X_CKFiterated(:,end)' + dx_CKFiterated(:,end)';
% 
% CKFtrace = zeros(1,length(t));
% CKFtrace_r = zeros(1,length(t));
% CKFtrace_v = zeros(1,length(t));
% CKFtraceIterated = zeros(1,length(t));
% for i = 1:length(t)
%     CKFtrace(i) = trace(P_CKF(1:6,1:6,i));
%     CKFtrace_r(i) = trace(P_CKF(1:3,1:3,i));
%     CKFtrace_v(i) = trace(P_CKF(4:6,4:6,i));
%     CKFtraceIterated(i) = trace(P_CKFiterated(1:6,1:6,i));
% end
% figure
% semilogy(t/60^2, CKFtrace)
% title("Trace of the Covariance of Spacecraft State")
% xlabel("time [hours]")
% ylabel("trace(P)  [km^2 + (km/s)^2]")
% figure
% semilogy(t/60^2, CKFtraceIterated)
% title("Trace of the Covariance of Spacecraft State (Iterated CKF)")
% xlabel("time [hours]")
% ylabel("trace(P)  [km^2 + (km/s)^2]")



