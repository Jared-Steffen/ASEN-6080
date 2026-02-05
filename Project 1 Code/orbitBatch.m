function [Xhat,P,y,yhat,batch_cnt] = orbitBatch(t,xbar0,Pbar0,Y,R,Xnom,constants,stations,tol)
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
    >tol: error tolerance to stop iterations of batch filter
Outputs:
    >Xhat: full predicted state
    >P: state history of estimated state covariance
    >y: pre-fit residuals for each batch
    >yhat: post-fit residuals for each batch
    >batch_cnt: number of batches it took to converge
%}

% Extract info
station_ids = stations.station_ids;

% Initialize % edit this!!!!
xbar = xbar0;
batch_cnt = 0;
err = 0.1;
n = length(xbar0);
Xhat(:,:,1) = Xnom;
nt = length(t);
P = zeros(n,n,nt);

% Iterative algorithm
while err > tol
    k = batch_cnt + 1;

    % Reset for iteration
    PbarR = chol(Pbar0);
    PbarR_inv = inv(PbarR);
    Pbar_inv = PbarR_inv*PbarR_inv';
    Lambda = Pbar_inv;
    N = Lambda*xbar;
    Phii = eye(n);
    Xbari = Xnom(1,:);

    % Integrate current initial condition
    state_stm = [Xbari';reshape(Phii,[],1)];
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
    [~,sstm] = ode45(@(t,x) odeSTM_J2_Drag(t,x,constants),t,state_stm,options);
    Xhat(:,:,k) = sstm(:,1:n);

    for i = 1:length(t)

        % Get this time step's state info
        Xbari(i,:) = sstm(i,1:n);
        Phii(:,:,i) = reshape(sstm(i,n+1:end),n,n);
    
        % Read next observation and determine station to pass to measurements
        M = Y(i,:);
        current_id = M(2);
        idx = find(current_id == station_ids);
        current_station = [current_id, Xbari(i,9+3*(idx-1)+(1:3))];

        % Extract station id and generate measurement
        [G,gs_state(:,i)] = genSingleMeasurement(t(i),stations,current_station,constants,Xbari(i,1:6));

        % Prefits and update
        y(:,i,k) = M(:,3:4)' - G(:,2:3)';
        Htilde = linearizedH(t(i),Xbari(i,1:6),[current_id;gs_state(:,i)],constants,station_ids);
        Hi = Htilde*Phii(:,:,i);
        Lambda = Lambda + Hi'/R*Hi;
        N = N + Hi'/R*y(:,i,k);
    end

    % Solve normal equations  
    LR = chol(Lambda);
    LR_inv = inv(LR);
    P0 = LR_inv*LR_inv';
    dxhat = P0*N;
    err = norm(dxhat);
    Xnom(1,:) = Xnom(1,:) + dxhat';
    xbar = xbar - dxhat;

    % Get covariance at each step and post fits
    for i = 1:length(t)
        P(:,:,i,k) = Phii(:,:,i)*P0*Phii(:,:,i)';
        M = Y(i,:);
        current_id = M(2);
        Htilde = linearizedH(t(i),Xbari(i,1:6),[current_id;gs_state(:,i)],constants,station_ids);
        Hi = Htilde*Phii(:,:,i);
        yhat(:,i,k) = y(:,i,k) - Hi*dxhat;
    end

    % Increase counter
    batch_cnt = batch_cnt + 1;

end

% Transpose outputs
y = pagetranspose(y);
yhat = pagetranspose(yhat);

end