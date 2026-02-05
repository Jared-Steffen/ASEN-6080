function [Xhat,P,y,yhat] = orbitLKF(t,xbar0,Pbar0,Y,R,Xnom,constants,stations,num_iterations)
%{
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Y: simulated actual measurement data for true trajectory
    >R: measurement noise covariance matrix
    >Xnom: nominal trajectory state vector
    >constants: struct that contains:
        -mu: gravitational parameter for central body
        -Rp: radius of central body
        -J2: J2 coefficient of central body
    >stations: struct containing information on stations:
        -Rs: positions of each station in the ecef frame
        -station_ids: stations ids correspoding to stations in Rs
        -theta0: initial spin angle of Earth in simulation in radians
        -el_mask: elevation mask for GS to S/C visibility in radians
    >num_iterations: number of times to iterate LKF algorithm

Outputs:
    >Xhat: full predicted state (dxhat + Xnom)
    >P: state history of estimated state covariance
    >y: pre-fit measurement residuals
    >yhat: post-fit residuals
%}

% Extract info
station_ids = stations.station_ids;

% Initialize
n = length(xbar0);

Phi = zeros(n, n, length(t));
P = zeros(n, n, length(t), num_iterations);
y = NaN(2, length(t), num_iterations);
yhat = NaN(2, length(t), num_iterations);
Xhat = zeros(length(t), n, num_iterations);
X0 = Xnom(1,:);

for j = 1:num_iterations
    
    % Reset for this iteration
    Pim1 = Pbar0;
    xhatim1 = xbar0;
    
    % Integrate current initial condition
    state_stm = [X0';reshape(eye(n),[],1)];
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
    [~,sstm] = ode45(@(t,x) odeSTM_J2_Drag(t,x,constants),t,state_stm,options);
    Xbar = sstm(:,1:n);


    % Iterative algorithm
    for i = 1:length(t)

        % Get this time step's state info
        Xbari(i,:) = Xbar(i,:);
        Phi(:,:,i) = reshape(sstm(i,n+1:end),n,n);
    
        % Read next observation and determine station to pass to measurements
        M = Y(i,:);
        current_id = M(2);
        idx = find(current_id == station_ids);
        current_station = [current_id, Xbari(i,9+3*(idx-1)+(1:3))];
    
        % Time update
        if i > 1
            Phii = Phi(:,:,i)*Phi(:,:,i-1)\eye(n);
        else
            Phii = Phi(:,:,i);
        end
        xbari = Phii*xhatim1;
        Pbari = Phii*Pim1*Phii';

        % Extract station id and generate measurement
        [G,gs_state(:,i)] = genSingleMeasurement(t(i),stations,current_station,constants,Xbari(i,1:6));

        yi = M(:,3:4)' - G(:,2:3)';
        Htilde = linearizedH(t(i),Xbari(i,1:6),[current_id;gs_state(:,i)],constants,station_ids);
        Ki = Pbari*Htilde'/(Htilde*Pbari*Htilde' + R);
    
        % Measurement correction
        y(:,i,j) = yi - Htilde*xbari;
        dxhat(:,i) = xbari + Ki*y(:,i,j);
        P(:,:,i,j) = (eye(n) - Ki*Htilde)*Pbari*(eye(n) - Ki*Htilde)' + Ki*R*Ki';
        yhat(:,i,j) = yi - Htilde*dxhat(:,i);
    
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

