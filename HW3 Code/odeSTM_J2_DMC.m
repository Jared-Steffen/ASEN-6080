function [stm_var_dot] = odeSTM_J2_DMC(~,stm_var,constants)
%{
NOTE: This considers the effects of J2 orbital pertubations
Inputs:
    >stm_state_vec: n^2 + n length vector of all STM entries and state vector
    >constants: struct of constants that must contain
        -mu: Earth's gravitational parameter mu
        -J2: Earth's J2 coefficient
        -RE: Earth's radius
Outputs:
    >stm_state_vec_dot: t.r.o.c of stm_state_vec
%}

% Extract constants
mu = constants.mu;
J2 = constants.J2;
RE = constants.RE;
tau = constants.tau;

% Reshape STM and extract necessary components
x = stm_var(1);
y = stm_var(2);
z = stm_var(3);
vx = stm_var(4);
vy = stm_var(5);
vz = stm_var(6);
wx = stm_var(7);
wy = stm_var(8);
wz = stm_var(9);
r = norm(stm_var(1:3));
phi = reshape(stm_var(10:end),6,6);

% Assign t.r.o.c variables
x_dot = vx;
y_dot = vy;
z_dot = vz;
vx_dot = (-mu*x)/r^3 + (1.5*mu*RE^2*J2)/r^5*(5*z^2/r^2-1)*x + wx;
vy_dot = (-mu*y)/r^3 + (1.5*mu*RE^2*J2)/r^5*(5*z^2/r^2-1)*y + wy;
vz_dot = (-mu*z)/r^3 + (1.5*mu*RE^2*J2)/r^5*(5*z^2/r^2-3)*z + wz;
wx_dot = -1/tau*wx;
wy_dot = -1/tau*wy;
wz_dot = -1/tau*wz;

% Initialize
A = zeros(6);

% Fill in dv/dv portion
A(1:3,4:6) = eye(3);

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

A(4:6,1:3) = da_mudr + da_J2dr;

phi_dot = A*phi;

state_vec_dot = [x_dot;y_dot;z_dot;vx_dot;vy_dot;vz_dot;wx_dot;wy_dot;wz_dot];

stm_var_dot = [state_vec_dot;reshape(phi_dot,[],1)];


end
