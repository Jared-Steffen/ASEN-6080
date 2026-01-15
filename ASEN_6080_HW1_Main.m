clc; clear; close all

%% Problem 1

% Read in data
file_p1 = fileread('HW1 Data\prob1c_solution.json');
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
file_p2 = fileread('HW1 Data\prob2b_solution.json');
data_p2 = jsondecode(file_p2);
% disp(data_p2)

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
mu = 398600;  % km^3/s^2
J2 = 1.08264e-3;

% Get initial r and v
[r0,v0] = oe2rv(mu,a,e,Omega,i,w,nu0);
var = [r0;v0;J2];
pert = [1 0 0 0 0.001 0 0]';
var_pert = var + pert;

% Orbital period
n = sqrt(mu/a^3);
T = (2*pi)/n;

% Time vector
tspan = [0 15*T];

% ode45 calls for truth data
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,state_unpert] = ode45(@(t,x) orbitEOM_J2(t,x,mu,Re,J2),tspan,var,options);
[~,state_pert] = ode45(@(t,x) orbitEOM_J2(t,x,mu,Re,J2),t,var_pert,options);

deltaxs{1} = state_pert - state_unpert;
times{1} = t;

% Reshape ICs for STM
% var_stm = [test_X0;reshape(test_Phi0,[],1)];
var_stm = [var_pert;reshape(eye(7),[],1)];

tspan = linspace(0,15*T,length(t));

% ode45 calls for STM
[t,stm_pert] = ode45(@(t,x) odeSTM_J2(t,x,mu,Re),tspan,var_stm,options);

% Extract STMs and reshape to be 7x7 and then propogate pertubations
for i = 1:length(t)
    stm_p2 = reshape(stm_pert(i,8:end),7,7);
    pertt(i,:) = stm_p2 * pert;
end

deltaxs{2} = pertt;
times{2} = t;

% Plots
plot_rv_state(times,deltaxs,{"NL Propagation","STM Propagation"},true)

r_labels = {'\deltax Position Difference [km]','\deltay Position Difference [km]','\deltaz Position Difference [km]'};
v_labels = {'\deltav_x Velocity Difference [km/s]','\deltav_y Velocity Difference [km/s]','\deltav_z Velocity Difference [km/s]'};

pert_diff = (state_pert - state_unpert) - pertt;
r = pert_diff(:,1:3);
v = pert_diff(:,4:6);

figure();
subplot(311)
hold on
plot(t/3600,r(:,1))
xlabel('Time [hr]')
ylabel(r_labels{1})
subplot(312)
hold on
plot(t/3600,r(:,2))
xlabel('Time [hr]')
ylabel(r_labels{2})
subplot(313)
hold on
plot(t/3600,r(:,3))
xlabel('Time [hr]')
ylabel(r_labels{3})
sgtitle('\deltar_{NL} - \deltar_{STM} Position Difference')

figure();
subplot(311)
hold on
plot(t/3600,v(:,1))
xlabel('Time [hr]')
ylabel(v_labels{1})
subplot(312)
hold on
plot(t/3600,v(:,2))
xlabel('Time [hr]')
ylabel(v_labels{2})
subplot(313)
hold on
plot(t/3600,v(:,3))
xlabel('Time [hr]')
ylabel(v_labels{3})
sgtitle('\deltav_{NL} - \deltav_{STM} Velocity Difference')

%% Problem 3
