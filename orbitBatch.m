function [Xhat,P,yhat,batch_cnt] = orbitBatch(t,xbar0,Pbar0,Y,R,Xnom,gs_state,mu,Rp,J2,measurement_params,tol)
%{
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Y: simulated actual measurement data for true trajectory
    >R: measurement noise covariance matrix
    >Xnom: nominal trajectory state vector
    >gs_state: state of the ground stations in eci frame at each time step
    >mu: gravitational parameter for central body
    >Rp: radius of central body
    >J2: J2 coefficient of central body
    >measurement_params: struct containing measurement station locations in
                         lla, initial angle of central body spin, angular
                         velocity of central body, elevation mask in
                         degrees
    >tol: error tolerance to stop iterations of batch filter
Outputs:
    >Xhat: full predicted state
    >P: state history of estimated state covariance
    >yhat: post-fit residuals for each batch
    >batch_cnt: number of batches it took to converge
%}

% Extract measurement params
stations_lla = measurement_params.sall_lla;
theta0 = measurement_params.theta0;
wE = measurement_params.wE;
el_mask = measurement_params.el_mask;

% Initialize
xbar = xbar0;
batch_cnt = 0;
err = 0.1;
n = length(xbar0);
Xhat(:,:,1) = Xnom;

% Iterative algorithm
while err > tol
    k = batch_cnt + 1;
    Lambda = inv(Pbar0);
    N = Lambda*xbar;
    Phii = eye(n);
    Xbari = Xnom(1,:);
    for i = 1:length(t)
    
        % Read next observation and expected observation
        M = Y{i};
    
        % Time update
        if i > 1
            current_tstep = [t(i-1) t(i)];
            state_stm = [Xbari';reshape(Phii,[],1)];
            options = odeset('RelTol',1e-10,'AbsTol',1e-10);
            [~,sstm] = ode45(@(t,x) odeSTM_J2_rv(t,x,mu,Rp,J2),current_tstep,state_stm,options);
            last_sstm = sstm(end,:);
            Xbari = last_sstm(1:n);
            Phii = reshape(last_sstm(n+1:end),n,n);
        end
        Xhat(i,:,k) = Xbari';
        Phii_store(:,:,i) = Phii;
    
        % Computes observation deviation
        if isempty(M)
            yhat(:,i,k) = NaN;
            continue % skips measurement update if no measurement available
        end
        [G,~] = genMeasurements(t(i),stations_lla,theta0,wE,el_mask,Xbari);
        Mn = G{1};
        if isempty(Mn)
            continue
        end
        yhat(:,i,k) = M(:,2:3)' - Mn(:,2:3)';
        Htilde = sc_range_ranger_Htilde(Xbari,gs_state(i,:));
        Hi = Htilde*Phii;
        Lambda = Lambda + Hi'/R*Hi;
        N = N + Hi'/R*yhat(:,i,k);
    end

    % Get covariance at each step
    for i = 1:length(t)
        P(:,:,i,k) = Phii_store(:,:,i)/Lambda*Phii_store(:,:,i)';
    end

    % Solve normal equations
    dxhat = Lambda\N;
    err = norm(dxhat);
    Xnom(1,:) = Xnom(1,:) + dxhat';
    xbar = xbar - dxhat;

    % Increase counter
    batch_cnt = batch_cnt + 1;
end

% Transpose outputs
yhat = pagetranspose(yhat);

end