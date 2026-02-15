function [measurements] = genMeasurements(t,stations,constants,sc_state)
%{
Inputs:
    >t: time index of measurement to be generates
    >stations: struct containing information on stations:
        -lla: positions of each station in lla coordinates
        -Rs: positions of each station in the ecef frame
        -station_ids: stations ids correspoding to stations in Rs
        -el_mask: elevation mask for GS to S/C visibility in radians
    >constants: struct of constants that must contain
        -wE: Earth's rotation rate
    >sc_state: Nx6 state vector of S/C through the simulation in inertial
               frame at time t
Outputs:
    >measurements: cell array of measurements at each time step in the
                    form of [t, station id, range, range rate]
                    with units [s, -, km, km/s]
%}

% Extract constants
wE = constants.wE;
el_mask = stations.el_mask;

% Extract station id and state in ecef
Rs = stations.Rs;

% Create ECI Earth angular speed vector
wE_eci = [0 0 wE]';

% Initialize
N = length(t);
k = 1;

% Iterate through each time step
for i = 1:N
    % Model spinning of Earth
    theta_t = wE*t(i) + constants.theta0;

    % Get S/C position and velocity for time step
    sc_pos_eci = sc_state(i,1:3)';
    sc_vel_eci = sc_state(i,4:6)';

    % Rotation matrix from ECEF to ECI
    R3 = [cos(theta_t) -sin(theta_t) 0;
          sin(theta_t) cos(theta_t) 0;
          0 0 1];
    
    % Calculate elevation and apply mask for taking measurements
    for j = 1:size(Rs,1)

        % Rotate stations to ECI frame, get velocity
        r_station_eci = R3 * Rs(j,:)';
        v_station_eci = cross(wE_eci,r_station_eci);
        
        % Rotate S/C position into ECEF
        sc_pos_ecef = R3'*sc_pos_eci;

        % Get LLA coords
        station_lla = stations.lla(j,:);

        % Rotation matrix from ECEF to Topocentric
        Rtop = [-sind(station_lla(2)), cosd(station_lla(2)), 0;
                -sind(station_lla(1))*cosd(station_lla(2)),...
                -sind(station_lla(1))*sind(station_lla(2)),...
                cosd(station_lla(1));
                cosd(station_lla(1))*cosd(station_lla(2)),...
                cosd(station_lla(1))*sind(station_lla(2)),...
                sind(station_lla(1))];
            
        % Rotate range into Topocentric
        rR_ecef = sc_pos_ecef - Rs(j,:)';
        rR_top = Rtop*rR_ecef;

        % Determine elevation
        el_angle = atan2(rR_top(3),norm(rR_top(1:2)));

        % Determine if measurements for station j are recorded
        if el_angle >= el_mask

            % Record measurements
            rho = norm(sc_pos_eci - r_station_eci);
            rho_dot = dot(sc_pos_eci-r_station_eci,sc_vel_eci-v_station_eci)/rho;
            measurements(k,:) = [t(i); stations.station_ids(j); rho; rho_dot];

            % Increase counter
            k = k+1;

        end
    end  
end


end

