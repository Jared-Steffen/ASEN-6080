function [var_dot] = orbitEOM_SunEarth3BP_SRP(t,var,constants)

% Extract constants
muE = constants.muEarth;
muS = constants.muSun;
P_Phi = constants.P_Phi;
AoM = constants.srpAMratio;

% Extract state variables
x = var(1);
y = var(2);
z = var(3);
vx = var(4);
vy = var(5);
vz = var(6);
CR = var(7);

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
rS_sc_hat = rS_sc_vec/norm(rS_sc_vec);

% Assign t.r.o.c variables
K = CR*AoM*P_Phi/1000;
x_dot = vx;
y_dot = vy;
z_dot = vz;
v_dot = -muE/r^3*r_vec + muS/rS_sc^3*rS_sc_vec - muS/RE_S^3*RS_E_vec - K*(constants.kmAU/rS_sc)^2.*rS_sc_hat;

var_dot = [x_dot;y_dot;z_dot;v_dot(:);0];

end

