clc; clear; close all

%% Problem 1

% Read in data
file = fileread('HW1 Data\prob1c_solution.json');
data = jsondecode(file);
disp(data)

% Constants
mu = data.inputs.state.mu;
J2 = data.inputs.state.J2;
J3 = data.inputs.state.J3;
r_vec = data.inputs.state.r;
Re = 6378; % km

% Outputs
A_ans = data.outputs.A_matrix.values;

% STM fxn call
stm_jacobian = STM_J2_J3(mu,J2,J3,Re,r_vec);

