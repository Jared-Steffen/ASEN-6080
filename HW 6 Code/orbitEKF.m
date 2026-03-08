function [Xhat,P,y,yhat] = orbitEKF(t,xbar0,Pbar0,Q,Qframe,Y,R,Xnom,constants,stations,LKFinit)
%{
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Q: process noise covariance matrix
    >Qframe: frame that Q was defined in (ECI or RIC)
    >Y: simulated actual measurement data for true trajectory
    >R: measurement noise covariance matrix
    >constants: struct that contains:
        -mu: Earth's gravitational parameter mu
        -J2: Earth's J2 coefficient
        -RE: Earth's radius
    >stations: struct containing information on stations:
        -Rs: positions of each station in the ecef frame
        -station_ids: stations ids correspoding to stations in Rs
        -theta0: initial spin angle of Earth in simulation in radians
        -el_mask: elevation mask for GS to S/C visibility in radians
    >LKFinit: number of measurements to process with the LKF to initialize
Outputs:
    >Xhat: full predicted state (dxhat + Xnom)
    >P: state history of estimated state covariance
    >y: pre-fit measurement residuals
    >yhat: post-fit residuals
%}

% Set ode tolerance
options = odeset('RelTol',1e-11,'AbsTol',1e-11);

% Extract info
station_ids = stations.station_ids;

% Initialize
Q0 = Q;
n = length(xbar0);
Phi = zeros(n,n,length(t));
Pim1 = Pbar0;
xhatim1 = xbar0;
X0 = Xnom(1,:);

if LKFinit > 0
    % Determine time step to stop LKF initialization at
    tLKF = t(1:LKFinit);
    
    % Integrate current initial condition
    Xref = [X0';reshape(eye(n),[],1)];
    [~,Xint] = ode45(@(t,x) odeSTM_J2(t,x,constants),tLKF,Xref,options);
    Xbar = Xint(:,1:n);
    
    % Iterative algorithm
    for i = 1:length(tLKF)
    
        % Get this time step's state info
        Xbari = Xbar(i,:);
        Phi(:,:,i) = reshape(Xint(i,n+1:end),n,n);
    
        % Read next observation and determine station to pass to measurements
        M = Y(i,:);
        current_id = M(2);
        idx = find(current_id == station_ids);
        current_station = [current_id, stations.Rs(idx,:)];
    
        % Time update
        if i > 1
            Phii = Phi(:,:,i)/Phi(:,:,i-1);
            deltat = t(i) - t(i-1);
        else
            Phii = Phi(:,:,i);
            deltat = 0;
        end
        
        % Rotate from RIC to ECI
        if ~isempty(Q)
            Q = Q0;
            if Qframe == "RIC"
                Rbar = Xbari(1:3)./norm(Xbari);
                Cbar = cross(Xbari(1:3),Xbari(4:6))./norm(cross(Xbari(1:3),Xbari(4:6)));
                Ibar = cross(Cbar,Rbar);
                R_ECI2RIC = [Rbar; Ibar; Cbar];
                Q = R_ECI2RIC'*Q*R_ECI2RIC;
            end
    
            % SNC implementation
            Gamma = deltat .* [deltat/2.*eye(3);
                                eye(3)];
    
            % Skip SNC if measurement gap is big enough
            if deltat > 10
                Pbari = Phii*Pim1*Phii';
            else
                Pbari = Phii*Pim1*Phii' + Gamma*Q*Gamma';
            end
        else
            Pbari = Phii*Pim1*Phii';
        end

        xbari = Phii*xhatim1;

        % Extract station id and generate measurement
        [G,gs_state] = genSingleMeasurement(t(i),current_station,constants,Xbari(1:6));

        y(:,i) = M(:,3:4)' - G(:,2:3)';
        Htilde = linearizedH(Xbari(1:6),gs_state);
        Ki = Pbari*Htilde'/(Htilde*Pbari*Htilde' + R);
    
        % Measurement correction
        dxhat(:,i) = xbari + Ki*(y(:,i)- Htilde*xbari);
        P(:,:,i) = (eye(n) - Ki*Htilde)*Pbari*(eye(n) - Ki*Htilde)' + Ki*R*Ki';
        yhat(:,i) = y(:,i) - Htilde*dxhat(:,i);
    
        % Move iteration forward
        xhatim1 = dxhat(:,i);
        Pim1 = P(:,:,i);
    end
    
    % Add deviations to nominal trajectory for LKF section and initialize EKF
    Xhat = Xnom(1:length(tLKF),:)' + dxhat;
end



% Iterative algorithm for EKF
for i = LKFinit+1:length(t)

    if i == 1
        Xbar(i,:) = X0;
        Phii  = eye(n);
        deltat = 0;
    else
        % Integrate reference trajectory and stm to next time step
        current_tstep = [t(i-1) t(i)];
        Xref = [Xhat(:,end);reshape(eye(n),[],1)];
    
        [~,Xint] = ode45(@(t,x) odeSTM_J2(t,x,constants),current_tstep,Xref,options);
        Xint_end = Xint(end,:);
        Xbar(i,:) = Xint_end(1:n);
        Phii = reshape(Xint_end(n+1:end),n,n);
        deltat = t(i) - t(i-1);
    end

    % Read next observation and determine station to pass to measurements
    M = Y(i,:);
    current_id = M(2);
    idx = find(current_id == station_ids);
    current_station = [current_id, stations.Rs(idx,:)];

    % Rotate from RIC to ECI
    if ~isempty(Q)
        Q = Q0;
        if Qframe == "RIC"
            Rbar = Xbari(1:3)./norm(Xbari);
            Cbar = cross(Xbari(1:3),Xbari(4:6))./norm(cross(Xbari(1:3),Xbari(4:6)));
            Ibar = cross(Cbar,Rbar);
            R_ECI2RIC = [Rbar; Ibar; Cbar];
            Q = R_ECI2RIC'*Q*R_ECI2RIC;
        end
    
        % SNC implementation
        Gamma = deltat .* [deltat/2.*eye(3);
                            eye(3)];
    
        % Skip SNC if measurement gap is big enough
        if deltat > 10
            Pbari = Phii*Pim1*Phii';
        else
            Pbari = Phii*Pim1*Phii' + Gamma*Q*Gamma';
        end
    else
        Pbari = Phii*Pim1*Phii';
    end

    % Extract station id and generate measurement
    [G,gs_state] = genSingleMeasurement(t(i),current_station,constants,Xbar(i,1:6));

    y(:,i) = M(:,3:4)' - G(:,2:3)';
    Htilde = linearizedH(Xbar(i,1:6),gs_state);
    Ki = Pbari*Htilde'/(Htilde*Pbari*Htilde' + R);

    % Measurement correction
    dXhat = Ki*y(:,i);
    Xhat(:,i) = Xbar(i,:)' + dXhat;
    yhat(:,i) = y(:,i) - Htilde*dXhat;
    P(:,:,i) = (eye(n) - Ki*Htilde)*Pbari*(eye(n) - Ki*Htilde)' + Ki*R*Ki';

    % Move iteration forward
    Pim1 = P(:,:,i);

end

% Transpose vectors
Xhat = Xhat';
yhat = yhat';
y = y';

end