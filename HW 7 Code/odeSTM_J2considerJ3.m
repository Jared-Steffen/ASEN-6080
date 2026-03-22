function [stm_var_dot] = odeSTM_J2considerJ3(~,stm_var,constants)
%{
NOTE: This integrates the effects of J2 orbital pertubations and considers
J3
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

% Reshape STM and extract necessary components
x = stm_var(1);
y = stm_var(2);
z = stm_var(3);
vx = stm_var(4);
vy = stm_var(5);
vz = stm_var(6);
r = norm(stm_var(1:3));
psi = reshape(stm_var(7:end),7,7);
% phi = psi(1:6,1:6);
% theta = psi(1:6,7);

% Assign t.r.o.c variables
x_dot = vx;
y_dot = vy;
z_dot = vz;
vx_dot = (-mu*x)/r^3 + (1.5*mu*RE^2*J2)/r^5*(5*z^2/r^2-1)*x;
vy_dot = (-mu*y)/r^3 + (1.5*mu*RE^2*J2)/r^5*(5*z^2/r^2-1)*y;
vz_dot = (-mu*z)/r^3 + (1.5*mu*RE^2*J2)/r^5*(5*z^2/r^2-3)*z;

% Initialize
A = zeros(6);
B = zeros(6,1);

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

da_dJ3 = (mu*RE^3)/(2*r^5).*[5*x/r*(7*z^3/r^3-3*z/r);
                             5*y/r*(7*z^3/r^3-3*z/r);
                             25*z^4/r^4-30*z^2/r^2+3];

% Full Jacobian
A(4:6,1:3) = da_mudr + da_J2dr;

% Consider Jacobian
B(4:6) = da_dJ3;

% Propagate STM
psi_dot = [A B; zeros(1,7)]*psi;

% Propagate consider STM
% theta_dot = A*theta+B;

% Create state vector
state_vec_dot = [x_dot;y_dot;z_dot;vx_dot;vy_dot;vz_dot];

% Append state vector, STM, and consider STM
stm_var_dot = [state_vec_dot;reshape(psi_dot,[],1)];


end
