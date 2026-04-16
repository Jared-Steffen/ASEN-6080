function [Htilde] = linearizedH(sc_state,gs_state)
%{
Inputs:
    >t: time index of measurement to be generates
    >sc_state: 6x1 S/C position and velocity state vector
    >gs_meas_state: 7x1 GS position [2:4] and velocity [5:7] state vector and
               id [1] for station making measurement
    >constants: struct of constants that must contain
        -wE: Earth's rotation rate
    >station_ids: the list of station ids, in a consistent order
Outputs:
    >Htilde: jacobian of partials of rho and rho_dot wrt S/C state
%}

% Break apart S/C state vector
n = length(sc_state);
x = sc_state(1);
y = sc_state(2);
z = sc_state(3);
vx = sc_state(4);
vy = sc_state(5);
vz = sc_state(6);

% Break apart GS state vector
xs = gs_state(1);
ys = gs_state(2);
zs = gs_state(3);
vxs = gs_state(4);
vys = gs_state(5);
vzs = gs_state(6);

dx = x-xs;
dy = y-ys;
dz = z-zs;
dvx = vx - vxs;
dvy = vy - vys;
dvz = vz - vzs;
rho = norm(sc_state(1:3)'-gs_state(1:3));
N = dx*dvx + dy*dvy + dz*dvz;

Htilde = [dx/rho, dy/rho, dz/rho, 0, 0, 0, ; ... % Row 1
          dvx/rho-(dx*N)/rho^3, dvy/rho-(dy*N)/rho^3, dvz/rho-(dz*N)/rho^3,... % Row 2, Columns 1-3
          dx/rho, dy/rho, dz/rho]; % Row 2, Columns 4-6

Htilde = [Htilde, zeros(2,n-6)];

end

