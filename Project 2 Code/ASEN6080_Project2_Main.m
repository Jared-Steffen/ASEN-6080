clc; clear; close all

%% Load in 2a/b truth data and traj data

% Load data
data_2a = readmatrix("Project2a_Obs.txt");
data_2b = readmatrix("Project2b_Obs.txt");
data_truth_traj = load("Project2_Prob2_truth_traj_50days.mat");


%% Set constants and station locations

% Constants
constants.epochJD = 2456296.25; % JD
constants.muSun = 132712440017.987; % km^3/s^2
constants.muEarth = 3.98600432896939e5; % km^3/s^2
constants.P_Phi = 1357/299792458; % W/m^2 / m/s
constants.srpAMratio = 0.01; % m^2/kg
constants.kmAU = 149597870.7; % km
constants.wE = 7.29211585275553e-5; % rad/s
constants.RE = 6378.1363; % km
constants.theta0 = 0;

% Station info
stations.station_ids = [34, 65, 13];
rDSS34 = latlon2ECEF(constants.RE+0.691750,-35.398333,148.981944);
rDSS65 = latlon2ECEF(constants.RE+0.83453,40.427222,355.749444); % SECOND ANGLE NEGATIVE ONLY FOR PART 2
rDSS13 = latlon2ECEF(constants.RE+1.07114904,35.247164,243.205);
stations.Rs = [rDSS34';rDSS65';rDSS13'];

%% Flatten Obs Data (won't have to rework my code)

[t2a, obsData2a] = flattenObsData(data_2a,stations.station_ids);
[t2b, obsData2b] = flattenObsData(data_2b,stations.station_ids);
station_id2a = obsData2a(:,2);
station_id2b = obsData2b(:,2);

%% Filter Initialization -- Part a
X0 = [-274096790.0 -92859240.0 -40199490.0 32.67 -8.94 -3.88 1.2]';
X0_true = data_truth_traj.Xt_50(1,1:7);
n = length(X0);

P0 = diag([100 100 100 0.1 0.1 0.1 0.1].^2);

R = diag([5e-3 0.5e-6].^2);
V = chol(R,'lower');

%% Verify Dynamics and STM

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
t_test = data_truth_traj.Tt_50;

[~,test_stateSTM] = ode45(@(t,x) odeSTM_SunEarth3BP_SRP(t,x,constants),t_test,[X0(:);reshape(eye(n),[],1)],options);

diff = test_stateSTM(end,:) - data_truth_traj.Xt_50(end,:);

%% Truth trajectory
[~,Xtrue] = ode45(@(t,x) orbitEOM_SunEarth3BP_SRP(t,x,constants),t2a,X0_true(:),options);

%% Known State Estimate

% Nominal trajectory
[~,X_STM_nom_2a] = ode45(@(t,x) odeSTM_SunEarth3BP_SRP(t,x,constants),t2a,[X0(:);reshape(eye(n),[],1)],options);
Xnom_2a = X_STM_nom_2a(:,1:n);

% DCO Config
DCO_config = 0;
current_time = t2a;
current_data = obsData2a;
if DCO_config == 0 % All data
    DCO_cutoff_t = current_time;
    DCO_cutoff_data = current_data;  
    DCO_cutoff_Xnom = Xnom_2a;
    DCO_cutoff_stations = station_id2a;
elseif DCO_config == 1 % 50 days
    DCO_idxs = find(current_time < 50*86400);
    DCO_cutoff_t = current_time(DCO_idxs);
    DCO_cutoff_data = current_data(DCO_idxs,:);
    DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
    DCO_cutoff_stations = station_id2a(DCO_idxs);
elseif DCO_config == 2 % 100 days
    DCO_idxs = find(current_time < 100*86400);
    DCO_cutoff_t = current_time(DCO_idxs);
    DCO_cutoff_data = current_data(DCO_idxs,:);
    DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
    DCO_cutoff_stations = station_id2a(DCO_idxs);
elseif DCO_config == 3 % 150 days
    DCO_idxs = find(current_time < 150*86400);
    DCO_cutoff_t = current_time(DCO_idxs);
    DCO_cutoff_data = current_data(DCO_idxs,:);
    DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
    DCO_cutoff_stations = station_id2a(DCO_idxs);
else
    error('Error: Invalid DCO')
end


% Run filter
sigma_SNC = 1e-9;
Q = diag(sigma_SNC^2*ones(1,3));
% LKFinit = 10;
% [Xhat_2a,P_2a,y_2a,yhat_2a] = orbitEKF(DCO_cutoff_t,zeros(7,1),P0,[],"ECI",DCO_cutoff_data,R,DCO_cutoff_Xnom,constants,stations,LKFinit);
% plot_filter_diagnostics(DCO_cutoff_t,Xhat_2a,P_2a,y_2a,yhat_2a,R,DCO_cutoff_stations,"EKF","")
% [Xhat_2a,Xhatl_2a,dxhat_2a,dxhatl_2a,P_2a,Pl_2a,y_2a,yhat_2a] = orbitLKF(DCO_cutoff_t,zeros(7,1),P0,Q,"ECI",DCO_cutoff_data,V,DCO_cutoff_Xnom,constants,stations,3);
% plot_filter_diagnosticsV2(DCO_cutoff_t,Xtrue,Xhat_2a(:,:,end),P_2a(:,:,:,end),y_2a(:,:,end),yhat_2a(:,:,end),R,DCO_cutoff_stations,"LKF Filtered","")
% plot_filter_diagnosticsV2(DCO_cutoff_t,Xtrue,Xhatl_2a(:,:,end),Pl_2a(:,:,:,end),y_2a(:,:,end),yhat_2a(:,:,end),R,DCO_cutoff_stations,"LKF Smoothed","")
% [Xhat_2a,Xhatl_2a,dxhat_2a,dxhatl_2a,P_2a,Pl_2a,y_2a,yhat_2a] = SRIF(DCO_cutoff_t,zeros(7,1),P0,Q,DCO_cutoff_data,V,DCO_cutoff_Xnom,constants,stations,3);
% plot_filter_diagnosticsV2(DCO_cutoff_t,Xtrue,Xhat_2a(:,:,end),P_2a(:,:,:,end),y_2a(:,:,end),yhat_2a(:,:,end),R,DCO_cutoff_stations,"SRIF Filtered","")
% plot_filter_diagnosticsV2(DCO_cutoff_t,Xtrue,Xhatl_2a(:,:,end),Pl_2a(:,:,:,end),y_2a(:,:,end),yhat_2a(:,:,end),R,DCO_cutoff_stations,"SRIF Smoothed","")
% [Xhat_2a,P_2a,y_2a,yhat_2a,batch_cnt] = orbitBatch(DCO_cutoff_t,zeros(7,1),P0,DCO_cutoff_data,R,DCO_cutoff_Xnom,constants,stations,1e-6,"");
% plot_filter_diagnosticsV2(DCO_cutoff_t,Xtrue,Xhat_2a(:,:,end),P_2a(:,:,:,end),y_2a(:,:,end),yhat_2a(:,:,end),R,DCO_cutoff_stations,"Batch","")

% %% B-plane: Known State Estimate
% 
% % DCO Config
% DCO_config = 1;
% current_time = t2a;
% current_data = obsData2a;
% if DCO_config == 0 % All data
%     DCO_cutoff_t = current_time;
%     DCO_cutoff_data = current_data;  
%     DCO_cutoff_Xnom = Xnom_2a;
%     DCO_cutoff_stations = station_id2a;
% elseif DCO_config == 1 % 50 days
%     DCO_idxs = find(current_time < 50*86400);
%     DCO_cutoff_t = current_time(DCO_idxs);
%     DCO_cutoff_data = current_data(DCO_idxs,:);
%     DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
%     DCO_cutoff_stations = station_id2a(DCO_idxs);
% elseif DCO_config == 2 % 100 days
%     DCO_idxs = find(current_time < 100*86400);
%     DCO_cutoff_t = current_time(DCO_idxs);
%     DCO_cutoff_data = current_data(DCO_idxs,:);
%     DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
%     DCO_cutoff_stations = station_id2a(DCO_idxs);
% elseif DCO_config == 3 % 150 days
%     DCO_idxs = find(current_time < 150*86400);
%     DCO_cutoff_t = current_time(DCO_idxs);
%     DCO_cutoff_data = current_data(DCO_idxs,:);
%     DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
%     DCO_cutoff_stations = station_id2a(DCO_idxs);
% else
%     error('Error: Invalid DCO')
% end
% 
% [Xhat_2a,Xhatl_2a,dxhat_2a,dxhatl_2a,P_2a,Pl_2a,y_2a,yhat_2a] = SRIF(DCO_cutoff_t,zeros(7,1),P0,Q,DCO_cutoff_data,V,DCO_cutoff_Xnom,constants,stations,3);
% [cov_BPlane1,BdotRhat1,BdotThat1] = BPlaneAnalysis(DCO_cutoff_t,3,Xhatl_2a(:,:,end),Pl_2a(:,:,:,end),constants);
% 
% % DCO Config
% DCO_config = 2;
% current_time = t2a;
% current_data = obsData2a;
% if DCO_config == 0 % All data
%     DCO_cutoff_t = current_time;
%     DCO_cutoff_data = current_data;  
%     DCO_cutoff_Xnom = Xnom_2a;
%     DCO_cutoff_stations = station_id2a;
% elseif DCO_config == 1 % 50 days
%     DCO_idxs = find(current_time < 50*86400);
%     DCO_cutoff_t = current_time(DCO_idxs);
%     DCO_cutoff_data = current_data(DCO_idxs,:);
%     DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
%     DCO_cutoff_stations = station_id2a(DCO_idxs);
% elseif DCO_config == 2 % 100 days
%     DCO_idxs = find(current_time < 100*86400);
%     DCO_cutoff_t = current_time(DCO_idxs);
%     DCO_cutoff_data = current_data(DCO_idxs,:);
%     DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
%     DCO_cutoff_stations = station_id2a(DCO_idxs);
% elseif DCO_config == 3 % 150 days
%     DCO_idxs = find(current_time < 150*86400);
%     DCO_cutoff_t = current_time(DCO_idxs);
%     DCO_cutoff_data = current_data(DCO_idxs,:);
%     DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
%     DCO_cutoff_stations = station_id2a(DCO_idxs);
% else
%     error('Error: Invalid DCO')
% end
% 
% [Xhat_2a,Xhatl_2a,dxhat_2a,dxhatl_2a,P_2a,Pl_2a,y_2a,yhat_2a] = SRIF(DCO_cutoff_t,zeros(7,1),P0,Q,DCO_cutoff_data,V,DCO_cutoff_Xnom,constants,stations,3);
% [cov_BPlane2,BdotRhat2,BdotThat2] = BPlaneAnalysis(DCO_cutoff_t,3,Xhatl_2a(:,:,end),Pl_2a(:,:,:,end),constants);
% 
% % DCO Config
% DCO_config = 3;
% current_time = t2a;
% current_data = obsData2a;
% if DCO_config == 0 % All data
%     DCO_cutoff_t = current_time;
%     DCO_cutoff_data = current_data;  
%     DCO_cutoff_Xnom = Xnom_2a;
%     DCO_cutoff_stations = station_id2a;
% elseif DCO_config == 1 % 50 days
%     DCO_idxs = find(current_time < 50*86400);
%     DCO_cutoff_t = current_time(DCO_idxs);
%     DCO_cutoff_data = current_data(DCO_idxs,:);
%     DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
%     DCO_cutoff_stations = station_id2a(DCO_idxs);
% elseif DCO_config == 2 % 100 days
%     DCO_idxs = find(current_time < 100*86400);
%     DCO_cutoff_t = current_time(DCO_idxs);
%     DCO_cutoff_data = current_data(DCO_idxs,:);
%     DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
%     DCO_cutoff_stations = station_id2a(DCO_idxs);
% elseif DCO_config == 3 % 150 days
%     DCO_idxs = find(current_time < 150*86400);
%     DCO_cutoff_t = current_time(DCO_idxs);
%     DCO_cutoff_data = current_data(DCO_idxs,:);
%     DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
%     DCO_cutoff_stations = station_id2a(DCO_idxs);
% else
%     error('Error: Invalid DCO')
% end
% 
% [Xhat_2a,Xhatl_2a,dxhat_2a,dxhatl_2a,P_2a,Pl_2a,y_2a,yhat_2a] = SRIF(DCO_cutoff_t,zeros(7,1),P0,Q,DCO_cutoff_data,V,DCO_cutoff_Xnom,constants,stations,3);
% [cov_BPlane3,BdotRhat3,BdotThat3] = BPlaneAnalysis(DCO_cutoff_t,3,Xhatl_2a(:,:,end),Pl_2a(:,:,:,end),constants);
% 
% % DCO Config
% DCO_config = 0;
% current_time = t2a;
% current_data = obsData2a;
% if DCO_config == 0 % All data
%     DCO_cutoff_t = current_time;
%     DCO_cutoff_data = current_data;  
%     DCO_cutoff_Xnom = Xnom_2a;
%     DCO_cutoff_stations = station_id2a;
% elseif DCO_config == 1 % 50 days
%     DCO_idxs = find(current_time < 50*86400);
%     DCO_cutoff_t = current_time(DCO_idxs);
%     DCO_cutoff_data = current_data(DCO_idxs,:);
%     DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
%     DCO_cutoff_stations = station_id2a(DCO_idxs);
% elseif DCO_config == 2 % 100 days
%     DCO_idxs = find(current_time < 100*86400);
%     DCO_cutoff_t = current_time(DCO_idxs);
%     DCO_cutoff_data = current_data(DCO_idxs,:);
%     DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
%     DCO_cutoff_stations = station_id2a(DCO_idxs);
% elseif DCO_config == 3 % 150 days
%     DCO_idxs = find(current_time < 150*86400);
%     DCO_cutoff_t = current_time(DCO_idxs);
%     DCO_cutoff_data = current_data(DCO_idxs,:);
%     DCO_cutoff_Xnom = Xnom_2a(DCO_idxs,:);
%     DCO_cutoff_stations = station_id2a(DCO_idxs);
% else
%     error('Error: Invalid DCO')
% end
% 
% [Xhat_2a,Xhatl_2a,dxhat_2a,dxhatl_2a,P_2a,Pl_2a,y_2a,yhat_2a] = SRIF(DCO_cutoff_t,zeros(7,1),P0,Q,DCO_cutoff_data,V,DCO_cutoff_Xnom,constants,stations,3);
% [cov_BPlane4,BdotRhat4,BdotThat4] = BPlaneAnalysis(DCO_cutoff_t,3,Xhatl_2a(:,:,end),Pl_2a(:,:,:,end),constants);
% 
% plot_Bplane(cat(3,cov_BPlane1,cov_BPlane2,cov_BPlane3,cov_BPlane4),[BdotRhat1,BdotRhat2,BdotRhat3,BdotRhat4],[BdotThat1,BdotThat2,BdotThat3,BdotThat4],[50,100,150,200])

%% Unknown Issues, No Truth

% Recommended ICS
X0 = [-274096770.76544 -92859266.4499061 -40199493.6677441 32.6704564599943 -8.93838913761049 -3.87881914050316 1.0];

[~,X_STM_nom_2b] = ode45(@(t,x) odeSTM_SunEarth3BP_SRP(t,x,constants),t2b,[X0(:);reshape(eye(n),[],1)],options);
Xnom_2b = X_STM_nom_2b(:,1:n);

% DCO Config
DCO_config = 3;
current_time = t2b;
current_data = obsData2b;
if DCO_config == 0 % All data
    DCO_cutoff_t = current_time;
    DCO_cutoff_data = current_data;  
    DCO_cutoff_Xnom = Xnom_2b;
    DCO_cutoff_stations = station_id2b;
elseif DCO_config == 1 % 50 days
    DCO_idxs = find(current_time < 50*86400);
    DCO_cutoff_t = current_time(DCO_idxs);
    DCO_cutoff_data = current_data(DCO_idxs,:);
    DCO_cutoff_Xnom = Xnom_2b(DCO_idxs,:);
    DCO_cutoff_stations = station_id2b(DCO_idxs);
elseif DCO_config == 2 % 100 days
    DCO_idxs = find(current_time < 100*86400);
    DCO_cutoff_t = current_time(DCO_idxs);
    DCO_cutoff_data = current_data(DCO_idxs,:);
    DCO_cutoff_Xnom = Xnom_2b(DCO_idxs,:);
    DCO_cutoff_stations = station_id2b(DCO_idxs);
elseif DCO_config == 3 % 150 days
    DCO_idxs = find(current_time < 150*86400);
    DCO_cutoff_t = current_time(DCO_idxs);
    DCO_cutoff_data = current_data(DCO_idxs,:);
    DCO_cutoff_Xnom = Xnom_2b(DCO_idxs,:);
    DCO_cutoff_stations = station_id2b(DCO_idxs);
else
    error('Error: Invalid DCO')
end

% Initial Results, no additional params
% [Xhat_2b,Xhatl_2b,dxhat_2b,dxhatl_2b,P_2b,Pl_2b,y_2b,yhat_2b] = SRIF(DCO_cutoff_t,zeros(7,1),P0,Q,DCO_cutoff_data,V,X0,constants,stations,3);
% plot_filter_diagnostics(DCO_cutoff_t,Xhat_2b(:,:,end),P_2b(:,:,:,end),y_2b(:,:,end),yhat_2b(:,:,end),R,DCO_cutoff_stations,"SRIF Filtered","")
% plot_filter_diagnostics(DCO_cutoff_t,Xhatl_2b(:,:,end),Pl_2b(:,:,:,end),y_2b(:,:,end),yhat_2b(:,:,end),R,DCO_cutoff_stations,"SRIF Smoothed","")

% Consider Analysis
consider_labels = {
    '\mu_{Earth}'
    '\mu_{Sun}'
    'Station 34_x'
    'Station 34_y'
    'Station 34_z'
    'Station 65_x'
    'Station 65_y'
    'Station 65_z'
    'Station 13_x'
    'Station 13_y'
    'Station 13_z'
};
c = zeros(11,1);
Pcc = diag([1e3 1e3 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2].^2);

[Xhat_2bcc,dxhat_2bcc,P_2bcc,Pxc_2bcc,Pxx_2bcc,S_2bcc,y_2bcc,yhat_2bcc,Psi] = SRIFConsiderCov(DCO_cutoff_t,zeros(7,1),c,P0,Pcc,DCO_cutoff_data,V,X0,constants,stations);
plot_consider_effects(DCO_cutoff_t,P_2bcc,Pxx_2bcc,S_2bcc,Pcc,consider_labels)

% Add in grav parameters, station positions
sigma_SNC = 1e-9;
Q = diag(sigma_SNC^2*ones(1,3));

X0_new = [X0 constants.muEarth constants.muSun rDSS34' rDSS65' rDSS13'];

locked_station = 999;
sigma_muE = 1e-8;
sigma_muS = 1e-8;
if locked_station == 34
    P0_new = blkdiag(P0, diag([sigma_muE sigma_muS 1e-8 1e-8 1e-8 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2].^2));
elseif locked_station == 65
    P0_new = blkdiag(P0, diag([sigma_muE sigma_muS 1e-2 1e-2 1e-2 1e-8 1e-8 1e-8 1e-2 1e-2 1e-2].^2));
elseif locked_station == 13
    P0_new = blkdiag(P0, diag([sigma_muE sigma_muS 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-8 1e-8 1e-8].^2));
elseif locked_station == 999
    P0_new = blkdiag(P0, diag([sigma_muE sigma_muS 1e-8 1e-8 1e-8 1e-2 1e-2 1e-2 1e-8 1e-8 1e-8].^2));
else
    error('Error: Locked station must be "34", "65", or "13"')
end

[Xhat_2b,Xhatl_2b,dxhat_2b,dxhatl_2b,P_2b,Pl_2b,y_2b,yhat_2b,yhat1_2b ] = SRIF2(DCO_cutoff_t,zeros(18,1),P0_new,Q,DCO_cutoff_data,V,X0_new,constants,stations,3);
plot_filter_diagnosticsV3(DCO_cutoff_t,Xhat_2b(:,:,end),P_2b(:,:,:,end),y_2b(:,:,end),yhat_2b(:,:,end),R,DCO_cutoff_stations,"SRIF Filtered","")
plot_filter_diagnosticsV3(DCO_cutoff_t,Xhatl_2b(:,:,end),Pl_2b(:,:,:,end),yhat_2b(:,:,end),yhat1_2b(:,:,end),R,DCO_cutoff_stations,"SRIF Smoothed","")
[cov_BPlane1,BdotRhat1,BdotThat1] = BPlaneAnalysisV2(DCO_cutoff_t,3,Xhatl_2b(:,:,end),Pl_2b(:,:,:,end),constants);

% DCO Config
DCO_config = 0;
current_time = t2b;
current_data = obsData2b;
if DCO_config == 0 % All data
    DCO_cutoff_t = current_time;
    DCO_cutoff_data = current_data;  
    DCO_cutoff_Xnom = Xnom_2b;
    DCO_cutoff_stations = station_id2b;
elseif DCO_config == 1 % 50 days
    DCO_idxs = find(current_time < 50*86400);
    DCO_cutoff_t = current_time(DCO_idxs);
    DCO_cutoff_data = current_data(DCO_idxs,:);
    DCO_cutoff_Xnom = Xnom_2b(DCO_idxs,:);
    DCO_cutoff_stations = station_id2b(DCO_idxs);
elseif DCO_config == 2 % 100 days
    DCO_idxs = find(current_time < 100*86400);
    DCO_cutoff_t = current_time(DCO_idxs);
    DCO_cutoff_data = current_data(DCO_idxs,:);
    DCO_cutoff_Xnom = Xnom_2b(DCO_idxs,:);
    DCO_cutoff_stations = station_id2b(DCO_idxs);
elseif DCO_config == 3 % 150 days
    DCO_idxs = find(current_time < 150*86400);
    DCO_cutoff_t = current_time(DCO_idxs);
    DCO_cutoff_data = current_data(DCO_idxs,:);
    DCO_cutoff_Xnom = Xnom_2b(DCO_idxs,:);
    DCO_cutoff_stations = station_id2b(DCO_idxs);
else
    error('Error: Invalid DCO')
end

[Xhat_2b,Xhatl_2b,dxhat_2b,dxhatl_2b,P_2b,Pl_2b,y_2b,yhat_2b,yhat1_2b ] = SRIF2(DCO_cutoff_t,zeros(18,1),P0_new,Q,DCO_cutoff_data,V,X0_new,constants,stations,3);
[cov_BPlane2,BdotRhat2,BdotThat2] = BPlaneAnalysisV2(DCO_cutoff_t,3,Xhatl_2b(:,:,end),Pl_2b(:,:,:,end),constants);


%% Maneuver Estimation
% X0_new = [X0_new 0 0 0];
% sigma_dV = 1e-4;
% P0_new = blkdiag(P0_new, diag([sigma_dV sigma_dV sigma_dV].^2));
% constants.tm = 18675000;
% 
% [Xhat_2b,Xhatl_2b,dxhat_2b,dxhatl_2b,P_2b,Pl_2b,y_2b,yhat_2b,yhat1_2b ] = SRIF3(DCO_cutoff_t,zeros(21,1),P0_new,Q,DCO_cutoff_data,V,X0_new,constants,stations,3);
% plot_filter_diagnosticsV3(DCO_cutoff_t,Xhat_2b(:,:,end),P_2b(:,:,:,end),y_2b(:,:,end),yhat_2b(:,:,end),R,DCO_cutoff_stations,"SRIF Filtered","")
% plot_filter_diagnosticsV3(DCO_cutoff_t,Xhatl_2b(:,:,end),Pl_2b(:,:,:,end),yhat_2b(:,:,end),yhat1_2b(:,:,end),R,DCO_cutoff_stations,"SRIF Smoothed","")

% [cov_BPlane3,BdotRhat3,BdotThat3] = BPlaneAnalysisV3(DCO_cutoff_t,3,Xhatl_2b(:,:,end),Pl_2b(:,:,:,end),constants);
% plot_Bplane(cat(3,cov_BPlane1,cov_BPlane2,cov_BPlane3),[BdotRhat1,BdotRhat2,BdotRhat3],[BdotThat1,BdotThat2,BdotThat3],{'OD1 150 Days','OD1 200 Days','OD2 200 Days'})

X0_new = [X0 rDSS65' 0 0 0];
sigma_dV = 1e-4;
P0_new = blkdiag(P0, diag([1e-2 1e-2 1e-2 sigma_dV sigma_dV sigma_dV].^2));
constants.tm = 18675000;

[Xhat_2b,Xhatl_2b,dxhat_2b,dxhatl_2b,P_2b,Pl_2b,y_2b,yhat_2b,yhat1_2b ] = SRIF4(DCO_cutoff_t,zeros(13,1),P0_new,Q,DCO_cutoff_data,V,X0_new,constants,stations,5);
plot_filter_diagnosticsV4(DCO_cutoff_t,Xhat_2b(:,:,end),P_2b(:,:,:,end),y_2b(:,:,end),yhat_2b(:,:,end),R,DCO_cutoff_stations,"SRIF Filtered","")
plot_filter_diagnosticsV4(DCO_cutoff_t,Xhatl_2b(:,:,end),Pl_2b(:,:,:,end),yhat_2b(:,:,end),yhat1_2b(:,:,end),R,DCO_cutoff_stations,"SRIF Smoothed","")
[cov_BPlane3,BdotRhat3,BdotThat3] = BPlaneAnalysisV3(DCO_cutoff_t,3,Xhatl_2b(:,:,end),Pl_2b(:,:,:,end),constants);
plot_Bplane(cat(3,cov_BPlane1,cov_BPlane2,cov_BPlane3),[BdotRhat1,BdotRhat2,BdotRhat3],[BdotThat1,BdotThat2,BdotThat3],{'OD1 150 Days','OD1 200 Days','OD2 200 Days'})

%%
Xhatl_2b0 = Xhatl_2b(1,:,end);
fprintf('--- Final Initial State ---\n');

for i = 1:length(Xhatl_2b0)
    fprintf('X(%d) = %.8f\n', i, Xhatl_2b0(i));
end

fprintf('--- First Initial State ---\n');

for i = 1:length(X0_new)
    fprintf('X(%d) = %.8f\n', i, X0_new(i));
end

fprintf('--- Difference in Initial States ---\n');

for i = 1:length(X0_new)
    fprintf('X(%d) = %.8f\n', i,Xhatl_2b0(i)-X0_new(i));
end