function [Htilde] = linearizedH2(t,sc_state,gs_meas_state,constants,station_ids)
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
x = sc_state(1);
y = sc_state(2);
z = sc_state(3);
vx = sc_state(4);
vy = sc_state(5);
vz = sc_state(6);

% Model spinning of Earth
theta_t = constants.wE*t;

% Rotation matrix from ECEF to ECI
R3 = [cos(theta_t) -sin(theta_t) 0;
      sin(theta_t) cos(theta_t) 0;
      0 0 1];

% Loop associate correct ground station
S = length(station_ids);
prhopRs = zeros(1,3*S);
prhopdotRs = zeros(1,3*S);
current_id = gs_meas_state(1);
for i = 1:S
    if current_id == station_ids(i)
        xs = gs_meas_state(2);
        ys = gs_meas_state(3);
        zs = gs_meas_state(4);
        vxs = gs_meas_state(5);
        vys = gs_meas_state(6);
        vzs = gs_meas_state(7);

        % Assign differnce vars
        dx = x-xs;
        dy = y-ys;
        dz = z-zs;
        dvx = vx - vxs;
        dvy = vy - vys;
        dvz = vz - vzs;
        rho = norm(sc_state(1:3)'-gs_meas_state(2:4));
        N = dx*dvx + dy*dvy + dz*dvz;

        % Fill appropriate partials
        prhopRs_eci = [-dx/rho -dy/rho -dz/rho];
        prhopdotRs_eci = [(-dvx-constants.wE*dy)/rho+(dx*N)/rho^3,...
                          (-dvy+constants.wE*dx)/rho+(dy*N)/rho^3,...
                           -dvz/rho+(dz*N)/rho^3];
        % Rotate to ecef
        prhopRs(3*i-3+(1:3)) = prhopRs_eci*R3;
        prhopdotRs(3*i-3+(1:3)) = prhopdotRs_eci*R3;
        break
    end
end

Htilde = [dx/rho, dy/rho, dz/rho, 0, 0, 0, 0, 0, 0, ... % Row 1, Columns 1-9
          prhopRs; % Row 1, Columns 10-end
          dvx/rho-(dx*N)/rho^3, dvy/rho-(dy*N)/rho^3, dvz/rho-(dz*N)/rho^3,... % Row 2, Columns 1-3
          dx/rho, dy/rho, dz/rho, 0, 0, 0,... % Row 2, Columns 4-9
          prhopdotRs]; % Row 2, Columns 10-end

end
