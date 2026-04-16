function [measurement] = genMeasurementUKF(t,current_station,constants,sigma_pts)
%{
Inputs:
    >t: time index of measurement to be generates
    >stations: struct containing information on stations:
        -Rs: positions of each station in the ecef frame
        -station_ids: stations ids correspoding to stations in Rs
        -el_mask: elevation mask for GS to S/C visibility in radians
    >current_station: station information (id [1] and ecef state [2:end])
                      for station requesting measurement from
    >constants: struct of constants that must contain
        -wE: Earth's rotation rate
    >sigma_pts: 6x13 matrix of sigma points to be propagated
Outputs:
    >measurements: cell array of measurements at each time step in the
                    form of [station id, range, range rate, elevation with
                    units [-, km, km/s, degrees]
%}

% Extract constants
wE = constants.wE;

% Extract station id and state in ecef
current_station_Rs = current_station(2:end)';

% Create ECI Earth angular speed vector
wE_eci = [0 0 wE]';

% Model spinning of Earth
theta_t = wE*t + constants.theta0;

% Get S/C position and velocity for time step
sigma_pts_pos_eci = sigma_pts(1:3,:);
sigma_pts_vel_eci = sigma_pts(4:6,:);

% Rotation matrix from ECEF to ECI
R3 = [cos(theta_t) -sin(theta_t) 0;
      sin(theta_t) cos(theta_t) 0;
      0 0 1];

% Rotate stations to ECI frame, get velocity
r_station_eci = R3 * current_station_Rs;
v_station_eci = cross(wE_eci,r_station_eci);

% Record measurement
rho_vec = sigma_pts_pos_eci - r_station_eci;
rho = vecnorm(rho_vec,2,1);
rho_dot = dot(sigma_pts_pos_eci-r_station_eci,sigma_pts_vel_eci-v_station_eci,1)./rho;
measurement = [rho; rho_dot];

end

