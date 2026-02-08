function [stm_var_dot] = odeSTM_J2_Drag(t,stm_var,constants)
%{
NOTE: This considers the effects of J2 and drag orbital pertubations
Inputs:
    >stm_var: n^2 + n length vector of all STM entries and state vector
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
    >stm_var_dot: t.r.o.c of stm_state_vec
%}

% State vector length
n = 18;

% Extract constants
RE = constants.RE;
wE = constants.wE;
rho0 = constants.rho0;
r0 = constants.r0;
H = constants.H;
A = constants.A;
m = constants.m;

% Extract state variables
x = stm_var(1);
y = stm_var(2);
z = stm_var(3);
vx = stm_var(4);
vy = stm_var(5);
vz = stm_var(6);
mu = stm_var(7);
J2 = stm_var(8);
CD = stm_var(9);
phi = reshape(stm_var(n+1:end),n,n);

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
pos_vel_dot = [x_dot;y_dot;z_dot;vx_dot;vy_dot;vz_dot];

% Initialize
F = zeros(n);

% Fill in dv/dv portion
F(1:3,4:6) = eye(3);

% Fill in grav coefficients and CD section
F(4:6,7) = -1/r^3 .* [x, y, z]' + (3*J2*RE^2)/(2*r^5) .* [x*(5*z^2/r^2-1), y*(5*z^2/r^2-1), z*(5*z^2/r^2-3)]'; % mu
F(4:6,8) = (3*mu*RE^2)/(2*r^5) .* [x*(5*z^2/r^2-1), y*(5*z^2/r^2-1), z*(5*z^2/r^2-3)]'; % J2
F(4:6,9) = (-rho_atm*A*vrel)/(2*m) .* [vx+wE*y, vy-wE*x, vz]'; % CD

% Acceleration partials
da_mudr = mu/r^5.*[-(r^2-3*x^2), 3*x*y, 3*x*z;
                   3*x*y, -(r^2-3*y^2), 3*y*z;
                   3*x*z, 3*y*z, -(r^2-3*z^2)];

da_J2dr = [(3*J2*mu*RE^2)/(2*r^5)*(5*z^2/r^2*(1-7*x^2/r^2)-(1-5*x^2/r^2)),...
           (15*J2*mu*RE^2*x*y)/(2*r^7)*(1-7*z^2/r^2),...
           (15*J2*mu*RE^2*x*z)/(2*r^7)*(3-7*z^2/r^2);
           (15*J2*mu*RE^2*x*y)/(2*r^7)*(1-7*z^2/r^2),...
           (3*J2*mu*RE^2)/(2*r^5)*(5*z^2/r^2*(1-7*y^2/r^2)-(1-5*y^2/r^2)),...
           (15*J2*mu*RE^2*y*z)/(2*r^7)*(3-7*z^2/r^2);
           (15*J2*mu*RE^2*x*z)/(2*r^7)*(3-7*z^2/r^2),...
           (15*J2*mu*RE^2*y*z)/(2*r^7)*(3-7*z^2/r^2),...
           (3*J2*mu*RE^2)/(2*r^5)*(5*z^2/r^2*(3-7*z^2/r^2)-3*(1-5*z^2/r^2))];

K = (CD*A)/(2*m);

vrelx = vx + wE*y;
vrely = vy - wE*x;
vrelz = vz;

drho_dx = -rho_atm*x/(H*r);
drho_dy = -rho_atm*y/(H*r);
drho_dz = -rho_atm*z/(H*r);

dV_dvx = vrelx / vrel;
dV_dvy = vrely / vrel;
dV_dvz = vrelz / vrel;

dV_dx = -wE*dV_dvy;
dV_dy =  wE*dV_dvx;


da_Ddr = -K.*[drho_dx*vrel*vrelx + rho_atm*(dV_dx*vrelx + vrel*0), ...
              drho_dy*vrel*vrelx + rho_atm*(dV_dy*vrelx + vrel*wE), ...
              drho_dz*vrel*vrelx;
              drho_dx*vrel*vrely + rho_atm*(dV_dx*vrely - vrel*wE), ...
              drho_dy*vrel*vrely + rho_atm*(dV_dy*vrely + vrel*0), ...
              drho_dz*vrel*vrely;
              drho_dx*vrel*vrelz + rho_atm*(dV_dx*vrelz), ...
              drho_dy*vrel*vrelz + rho_atm*(dV_dy*vrelz), ...
              drho_dz*vrel*vrelz];

da_dv = -K*rho_atm.*[dV_dvx*vrelx + vrel*1, dV_dvy*vrelx + vrel*0, dV_dvz*vrelx + vrel*0;
                     dV_dvx*vrely + vrel*0,   dV_dvy*vrely + vrel*1,   dV_dvz*vrely + vrel*0; ...
                     dV_dvx*vrelz + vrel*0,   dV_dvy*vrelz + vrel*0,   dV_dvz*vrelz + vrel*1];

F(4:6,1:3) = da_mudr + da_J2dr + da_Ddr;
F(4:6,4:6) = da_dv;

% Integrate STM
phi_dot = F*phi;

state_vec_dot = [pos_vel_dot;zeros(n-length(pos_vel_dot),1)];

stm_var_dot = [state_vec_dot;reshape(phi_dot,[],1)];

end
