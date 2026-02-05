function [var_dot] = orbitEOM_J2_Drag(~,var,constants)
%{
Inputs:
    >var: 18x1 S/C position and velocity state vector + constants
    >constants: struct of constants that must contain
        -mu: Earth's gravitational parameter mu
        -J2: Earth's J2 coefficient
        -RE: Earth's radius
        -wE: Earth's rotation rate
        -CD: S/C drag coefficient
        -A: S/C cross sectional area (assuming spherical)
        -m: S/C mass
        -H, r0, rho0: atmospheric drag model constants
Outputs:
    >var_dot: t.r.o.c of state vector
%}

% Extract constants
mu = constants.mu;
J2 = constants.J2;
RE = constants.RE;
wE = constants.wE;
rho0 = constants.rho0;
r0 = constants.r0;
H = constants.H;
A = constants.A;
m = constants.m;
CD = constants.CD;

% Extract state variables
x = var(1);
y = var(2);
z = var(3);
vx = var(4);
vy = var(5);
vz = var(6);
N = length(var);

% Calculate radius and velocity
r = sqrt(x^2+y^2+z^2);
vrel = sqrt((vx+wE*y)^2+(vy-wE*x)^2+vz^2);

% Calculate air density rho
rho_atm = rho0*exp(-(r-r0)/H);

% J2 and Drag Acceleration Leading Terms
J2_term = (1.5*mu*RE^2*J2)/r^5;
D_term = (rho_atm*CD*A*vrel)/(2*m);

% Assign t.r.o.c variables
x_dot = vx;
y_dot = vy;
z_dot = vz;
vx_dot = (-mu*x)/r^3 + J2_term*(5*z^2/r^2-1)*x - D_term*(vx+wE*y);
vy_dot = (-mu*y)/r^3 + J2_term*(5*z^2/r^2-1)*y - D_term*(vy-wE*x);
vz_dot = (-mu*z)/r^3 + J2_term*(5*z^2/r^2-3)*z - D_term*vz;

% Final state derivative
var_dot = [x_dot;y_dot;z_dot;vx_dot;vy_dot;vz_dot;zeros(N-6,1)];

end