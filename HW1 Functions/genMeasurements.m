function [measurements,gs_state] = genMeasurements(t,stations_lla,theta0,wE,el_mask,sc_state)
%{
Inputs:
    >t: time vector of simulation
    >stations_lla: matrix where each col is the lla coordinates of each
                     station in [deg, deg, m]
    >theta0: initial spin angle of Earth in simulation in radians
    >wE: rotation rate of Earth in radians/s
    >el_mask: elevation mask for GS to S/C visibility in radians
    >sc_state: state vector of S/C through the simulation in inertial
                frame
Outputs:
    >measurements: cell array of measurements at each time step in the
                    form of [station id, range, range rate, elevation with
                    units [-, km, km/s, degrees]
    >gs_state: state of the ground station at every time step
%}

% Convert lla to ECEF
stations_ecef = lla2ecef(stations_lla')./1000;
stations_ecef = stations_ecef';

% Create ECI Earth angular speed vector
wE_eci = [0 0 wE]';

% Initialize
N = length(t);
gs_state = zeros(6, N);
measurements = cell(N,1);

% Iterate through each time step
for i = 1:N
    % Model spinning of Earth
    theta_t = wE*t(i) + theta0;

    % Get S/C position and velocity for time step
    sc_pos_eci = sc_state(i,1:3)';
    sc_vel_eci = sc_state(i,4:6)';

    % Rotation matrix from ECEF to ECI
    R3 = [cos(theta_t) -sin(theta_t) 0;
          sin(theta_t) cos(theta_t) 0;
          0 0 1];
    
    % Calculate elevation and apply mask for taking measurements
    station_measurements = [];
    for j = 1:size(stations_ecef,2)

        % Rotate stations to ECI frame, get velocity
        r_station_eci = R3 * stations_ecef(:,j);
        v_station_eci = cross(wE_eci,r_station_eci);
        
        % Rotate S/C position into ECEF
        sc_pos_ecef = R3'*sc_pos_eci;

        % Rotation matrix from ECEF to Topocentric
        Rtop = [-sind(stations_lla(2,j)), cosd(stations_lla(2,j)), 0;
                -sind(stations_lla(1,j))*cosd(stations_lla(2,j)),...
                -sind(stations_lla(1,j))*sind(stations_lla(2,j)),...
                cosd(stations_lla(1,j));
                cosd(stations_lla(1,j))*cosd(stations_lla(2,j)),...
                cosd(stations_lla(1,j))*sind(stations_lla(2,j)),...
                sind(stations_lla(1,j))];
            
        % Rotate range into Topocentric
        rR_ecef = sc_pos_ecef - stations_ecef(:,j);
        rR_top = Rtop*rR_ecef;

        % Determine elevation
        el_angle = atan2(rR_top(3),norm(rR_top(1:2)));

        % Determine if measurements for station j are recorded
        if el_angle >= el_mask

            % Record measurements
            rho = norm(sc_pos_eci - r_station_eci);
            rho_dot = dot(sc_pos_eci-r_station_eci,sc_vel_eci-v_station_eci)/rho;
            station_measurements = [station_measurements;
                                    j, rho, rho_dot, rad2deg(el_angle)];
            gs_state(:,i) = [r_station_eci;v_station_eci];

        else
            continue
        end
    end
    if size(gs_state, 2) < i
        gs_state(:,i) = [zeros(6,1)];
    end

    % Save measurements for time step
    measurements{i} = station_measurements;
   
end

% Transpose GS state vector
gs_state = gs_state';

end

