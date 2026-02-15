function [Xhat,P,y,yhat] = orbitLKF(t,xbar0,Pbar0,Q,Qframe,Y,R,Xnom,constants,stations,num_iterations)
%{
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Q: process noise covariance matrix
    >Qframe: frame that Q was defined in (ECI or RIC)
    >Y: simulated actual measurement data for true trajectory
    >R: measurement noise covariance matrix
    >Xnom: nominal trajectory state vector
    >constants: struct that contains:
        -mu: Earth's gravitational parameter mu
        -J2: Earth's J2 coefficient
        -RE: Earth's radius
        -wE: Earth's rotation rate
    >stations: struct containing information on stations:
        -Rs: positions of each station in the ecef frame
        -station_ids: stations ids correspoding to stations in Rs
        -theta0: initial spin angle of Earth in simulation in radians
    >num_iterations: number of times to iterate LKF algorithm

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
n = length(xbar0);

Phi = zeros(n,n,length(t));
P = zeros(n,n,length(t),num_iterations);
y = NaN(2, length(t),num_iterations);
yhat = NaN(2, length(t),num_iterations);
Xhat = zeros(length(t),n,num_iterations);
X0 = Xnom(1,1:n);

for j = 1:num_iterations
    
    % Reset for this iteration
    Pim1 = Pbar0;
    xhatim1 = xbar0;

    % Integrate current initial condition
    [~,Xint] = ode45(@(t,x) odeSTM_J2(t,x,constants),t,[X0';reshape(eye(n),[],1)],options);
    Xbar = Xint(:,1:n);

    % Iterative algorithm
    for i = 1:length(t)

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

        xbari = Phii*xhatim1;
        Pbari = Phii*Pim1*Phii' + Gamma*Q*Gamma';

        % Extract station id and generate measurement
        [G,gs_state] = genSingleMeasurement(t(i),current_station,constants,Xbari(1:6));

        y(:,i,j) = M(:,3:4)' - G(:,2:3)';
        Htilde = linearizedH(Xbari(1:6),gs_state);
        Ki = Pbari*Htilde'/(Htilde*Pbari*Htilde' + R);
    
        % Measurement correction
        dxhat(:,i) = xbari + Ki*(y(:,i,j)- Htilde*xbari);
        P(:,:,i,j) = (eye(n) - Ki*Htilde)*Pbari*(eye(n) - Ki*Htilde)' + Ki*R*Ki';
        yhat(:,i,j) = y(:,i,j) - Htilde*dxhat(:,i);
    
        % Move iteration forward
        xhatim1 = dxhat(:,i);
        Pim1 = P(:,:,i,j);
        
    end

    % Add deviations to nominal trajectory
    Xhat(:,:,j) = Xbar + dxhat';

    % Backwards solve for new xhat0
    x_hat0 = Phi(:,:,end)\dxhat(:,end);
    X0 = X0 + x_hat0';
    xbar0 = xbar0 - x_hat0;
end


% Trasponse residuals
y = pagetranspose(y);
yhat = pagetranspose(yhat);

end

