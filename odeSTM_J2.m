function [stm_state_vec_dot] = odeSTM_J2(~,stm_state_vec,mu,R)
%{
NOTE: This considers the effects of J2 orbital pertubations
Inputs:
    >stm_state_vec: n^2 + n length vector of all STM entries and state vector
    >mu: gravitational parameter of central body (km^3/s^2)
    >R: radius of central body
    >J2: J2 coefficient of central body
Outputs:
    >stm_state_vec_dot: t.r.o.c of stm_state_vec
%}

% Reshape STM and extract necessary components
x = stm_state_vec(1);
y = stm_state_vec(2);
z = stm_state_vec(3);
vx = stm_state_vec(4);
vy = stm_state_vec(5);
vz = stm_state_vec(6);
J2 = stm_state_vec(7);
r = norm(stm_state_vec(1:3));
phi = reshape(stm_state_vec(8:end),7,7);

% Assign t.r.o.c variables
x_dot = vx;
y_dot = vy;
z_dot = vz;
u_dot = (-mu*x)/r^3 + (1.5*mu*R^2*J2)/r^5*(5*z^2/r^2-1)*x;
v_dot = (-mu*y)/r^3 + (1.5*mu*R^2*J2)/r^5*(5*z^2/r^2-1)*y;
w_dot = (-mu*z)/r^3 + (1.5*mu*R^2*J2)/r^5*(5*z^2/r^2-3)*z;

% Initialize
A = zeros(7);

% Fill in dv/dv portion
A(1:3,4:6) = eye(3);

% Fill in grav coefficients section
A(4:6,7) = (3*mu*R^2)/(2*r^5) .* [x*(5*z^2/r^2-1), y*(5*z^2/r^2-1), z*(5*z^2/r^2-3)]'; % J2

% Acceleration partials

da_mudr = mu/r^5.*[-(r^2-3*x^2), 3*x*y, 3*x*z;
                   3*x*y, -(r^2-3*y^2), 3*y*z;
                   3*x*z, 3*y*z, -(r^2-3*z^2)];

da_J2dr = [(3*J2*mu*R^2)/(2*r^5)*(5*z^2/r^2*(1-7*x^2/r^2)-(1-5*x^2/r^2)),...
           (15*J2*mu*R^2*x*y)/(2*r^7)*(1-7*z^2/r^2),...
           (15*J2*mu*R^2*x*z)/(2*r^7)*(3-7*z^2/r^2);
           (15*J2*mu*R^2*x*y)/(2*r^7)*(1-7*z^2/r^2),...
           (3*J2*mu*R^2)/(2*r^5)*(5*z^2/r^2*(1-7*y^2/r^2)-(1-5*y^2/r^2)),...
           (15*J2*mu*R^2*y*z)/(2*r^7)*(3-7*z^2/r^2);
           (15*J2*mu*R^2*x*z)/(2*r^7)*(3-7*z^2/r^2),...
           (15*J2*mu*R^2*y*z)/(2*r^7)*(3-7*z^2/r^2),...
           (3*J2*mu*R^2)/(2*r^5)*(5*z^2/r^2*(3-7*z^2/r^2)-3*(1-5*z^2/r^2))];

A(4:6,1:3) = da_mudr + da_J2dr;

phi_dot = A*phi;

state_vec_dot = [x_dot;y_dot;z_dot;u_dot;v_dot;w_dot;0];

stm_state_vec_dot = [state_vec_dot;reshape(phi_dot,[],1)];


end
