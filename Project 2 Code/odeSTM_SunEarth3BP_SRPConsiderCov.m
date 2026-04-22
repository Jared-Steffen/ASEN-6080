function [stm_var_dot] = odeSTM_SunEarth3BP_SRPConsiderCov(t,stm_var,constants)
%{
NOTE: This considers the effects of Sun-Earth 3BP and SRP pertubations
Inputs:
    >stm_var: n^2 + n length vector of all STM entries and state vector
    >constants: struct of constants that must contain
        -mu: Earth's gravitational parameter mu
        -RE: Earth's radius
        -wE: Earth's rotation rate
Outputs:
    >stm_var_dot: t.r.o.c of stm_state_vec
%}

% State vector length
n = 7;

% Extract constants
muE = constants.muEarth;
muS = constants.muSun;
P_Phi = constants.P_Phi;
AoM = constants.srpAMratio;

% Extract state variables
x = stm_var(1);
y = stm_var(2);
z = stm_var(3);
vx = stm_var(4);
vy = stm_var(5);
vz = stm_var(6);
CR = stm_var(7);
Psi = reshape(stm_var(n+1:end),18,18);

% Calculate radius vectors and get components
r_vec = [x;y;z];
r = sqrt(x^2+y^2+z^2);
[RE_S_vec, ~, ~] = Ephem(constants.epochJD + t/86400,3,'EME2000');
RS_E_vec = -RE_S_vec(:);
RE_S = norm(RS_E_vec);
rS_sc_vec = RS_E_vec - r_vec;
rS_sc = norm(rS_sc_vec);
xS_sc = rS_sc_vec(1);
yS_sc = rS_sc_vec(2);
zS_sc = rS_sc_vec(3);
rS_sc_hat = rS_sc_vec./norm(rS_sc_vec);

% Assign t.r.o.c variables
K = CR*AoM*P_Phi/1000;
x_dot = vx;
y_dot = vy;
z_dot = vz;
v_dot = -muE/r^3.*r_vec + muS/rS_sc^3.*rS_sc_vec - muS/RE_S^3.*RS_E_vec - K*(constants.kmAU/rS_sc)^2.*rS_sc_hat;

pos_vel_dot = [x_dot;y_dot;z_dot;v_dot];

% Initialize
F = zeros(n);
B = zeros(n,11);

% Fill in dv/dv portion
F(1:3,4:6) = eye(3);

% Acceleration partials

da_mudr = muE/r^5.*[-(r^2-3*x^2), 3*x*y, 3*x*z;
                    3*x*y, -(r^2-3*y^2), 3*y*z;
                    3*x*z, 3*y*z, -(r^2-3*z^2)];

da_3Bdr = muS/rS_sc^5.*[-(rS_sc^2-3*xS_sc^2), 3*xS_sc*yS_sc, 3*xS_sc*zS_sc;
                        3*xS_sc*yS_sc, -(rS_sc^2-3*yS_sc^2), 3*yS_sc*zS_sc;
                        3*xS_sc*zS_sc, 3*yS_sc*zS_sc, -(rS_sc^2-3*zS_sc^2)];

da_SRPdr = (K/rS_sc^5)*constants.kmAU^2.*[rS_sc^2-3*xS_sc^2, -3*xS_sc*yS_sc, -3*xS_sc*zS_sc;
                                          -3*xS_sc*yS_sc, rS_sc^2-3*yS_sc^2, -3*yS_sc*zS_sc;
                                          -3*xS_sc*zS_sc, -3*yS_sc*zS_sc, rS_sc^2-3*zS_sc^2];

F(4:6,1:3) = da_mudr + da_3Bdr + da_SRPdr;

% CR Partial
da_dCR = -(AoM*P_Phi/1000)*(constants.kmAU/rS_sc)^2.*rS_sc_hat;
F(4:6,7) = da_dCR;

% Consider Partials
da_dmuE = -1/r^3.*r_vec;
da_dmuS = 1/rS_sc^3.*rS_sc_vec - 1/RE_S^3.*RS_E_vec;
B(4:6,1:2) = [da_dmuE,da_dmuS];

Psi_dot = [F B; zeros(11,18)]*Psi;

stm_var_dot = [pos_vel_dot;0;reshape(Psi_dot,[],1)];
end

