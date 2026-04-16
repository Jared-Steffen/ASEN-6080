function [Xhat,Xhatl,dxhat,dxhatl,P,Pl,y,yhat] = orbitLKF(t,xbar0,Pbar0,Q,Qframe,Y,R,Xnom,constants,stations,num_iterations)
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
    >Xhat: full predicted state (dxhat + Xnom)/(used in smoother) 
    >dxhat: state estimate correction (used in smoother)
    >P: state history of estimated state covariance (used in smoother)
    >Pbari: time update covariances (used in smoother)
    >y: pre-fit measurement residuals
    >yhat: post-fit residuals
    >Phii: STM from (t_k-1, tk) for each time step (used in smoother)
%}

% Set ode tolerance
options = odeset('RelTol',1e-11,'AbsTol',1e-11);

% Extract info
station_ids = stations.station_ids;

% Initialize
Q0 = Q;
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
    [~,Xint] = ode45(@(t,x) odeSTM_SunEarth3BP_SRP(t,x,constants),t,[X0';reshape(eye(n),[],1)],options);
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
            Phii(:,:,i,j) = Phi(:,:,i)/Phi(:,:,i-1);
            deltat = t(i) - t(i-1);
        else
            Phii(:,:,i,j) = Phi(:,:,i);
            deltat = 0;
        end
        
        % If adding process noise, then implement SNC
        if ~isempty(Q)
            Q = Q0;
            % Rotate from RIC to ECI
            if Qframe == "RIC"
                Rbar = Xbari(1:3)./norm(Xbari);
                Cbar = cross(Xbari(1:3),Xbari(4:6))./norm(cross(Xbari(1:3),Xbari(4:6)));
                Ibar = cross(Cbar,Rbar);
                R_ECI2RIC = [Rbar; Ibar; Cbar];
                Q = R_ECI2RIC'*Q*R_ECI2RIC;
            end
    
            % SNC implementation
            Gamma = zeros(n,3);
            Gamma(1:6,:) = deltat .* [deltat/2.*eye(3);
                                       eye(3)];
    
            % Skip SNC if measurement gap is big enough
            if deltat > 60
                Pbar(:,:,i,j) = Phii(:,:,i,j)*Pim1*Phii(:,:,i,j)';
            else
                Pbar(:,:,i,j) = Phii(:,:,i,j)*Pim1*Phii(:,:,i,j)' + Gamma*Q*Gamma';
            end
        else
            Pbar(:,:,i,j) = Phii(:,:,i,j)*Pim1*Phii(:,:,i,j)';
        end

        xbari = Phii(:,:,i,j)*xhatim1(:);

        % Extract station id and generate measurement
        [G,gs_state] = genSingleMeasurement(t(i),current_station,constants,Xbari(1:6));

        y(:,i,j) = M(:,3:4)' - G(:,2:3)';
        Htilde = linearizedH(Xbari,gs_state);
        Ki = Pbar(:,:,i,j)*Htilde'/(Htilde*Pbar(:,:,i,j)*Htilde' + R);
    
        % Measurement correction
        dxhat(:,i,j) = xbari + Ki*(y(:,i,j)- Htilde*xbari);
        P(:,:,i,j) = (eye(n) - Ki*Htilde)*Pbar(:,:,i,j)*(eye(n) - Ki*Htilde)' + Ki*R*Ki';
        yhat(:,i,j) = y(:,i,j) - Htilde*dxhat(:,i,j);
    
        % Move iteration forward
        xhatim1 = dxhat(:,i,j);
        Pim1 = P(:,:,i,j);
        
    end

    % Add deviations to nominal trajectory
    Xhat(:,:,j) = Xbar + dxhat(:,:,j)';
    
    [Xhatl(:,:,j),dxhatl(:,:,j),Pl(:,:,:,j)] = orbitLKFSmoother(t,dxhat(:,:,j),Xbar,P(:,:,:,j),Pbar(:,:,:,j),Phii(:,:,:,j));

    % Backwards solve for new xhat0
    x_hat0 = Phi(:,:,end)\dxhatl(:,end,j);
    X0 = X0 + x_hat0';
    xbar0 = xbar0 - x_hat0;
end



% Trasponse residuals
y = pagetranspose(y);
yhat = pagetranspose(yhat);

end

