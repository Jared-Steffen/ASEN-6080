clc; clear; close all

%% Setup

% HW J3 Data
data = load("simulation_dataJ2J3_test.mat");

% Split up data
t = data.t;
Xnom = data.Xnom;
noisy_measurements = data.noisy_measurements;
station_id = noisy_measurements(:,2);
R = data.R;
xbar0 = data.xbar0;
Pbar0 = data.Pbar0;
X0 = Xnom(1,:) + xbar0';
constants = data.constants;
stations = data.stations;

% Givens
Pbar0 = 10.* eye(6);
c = constants.J3;
Pcc0 = constants.J3^2;

%% Consider Analysis

% Consider LKF Filter
[Xhat,Xchat,dxhat,dxchat,P,Pxx,Pc,Psi,y,yhat] = orbitLKF_consider(t,0.*xbar0,c,Pbar0,Pcc0,noisy_measurements,R,Xnom,constants,stations);

% Plots
plot_consider_analysis(t,Xnom,Xhat,P,Pxx)

% Map results back to epoch
Psi_tf = Psi(:,:,end);
new_dx0 = Psi_tf(1:6,1:6)\dxhat(:,end);
new_Pc0 = Psi_tf\Pc(:,:,end)/Psi_tf';

% Integrate forward
options = odeset('RelTol',1e-11,'AbsTol',1e-11);
[~,Xint] = ode45(@(t,x) odeSTM_J2considerJ3(t,x,constants),t,[Xnom(1,:)'+new_dx0(1:6);reshape(eye(7),[],1)],options);
Xbar = Xint(:,1:6);

for i = 1:length(t)

    % Get this time step's state info
    Xbari = Xbar(i,:);
    Psi_new = reshape(Xint(i,7:end),7,7);
    Pc_new(:,:,i) = Psi_new*new_Pc0*Psi_new';

end
Pxx = Pc_new(1:6,1:6,:);

% Plots
plot_consider_analysis(t,Xnom,Xbar,[],Pxx)