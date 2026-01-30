function [Xhat,P,y,yhat,batch_cnt] = orbitBatch(t,xbar0,Pbar0,Y,R,Xnom,gs_state,mu,Rp,J2,measurement_params,tol)
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
    >y: pre-fit residuals for each batch
    >yhat: post-fit residuals for each batch
    >batch_cnt: number of batches it took to converge
%}

% Extract measurement params
stations_lla = measurement_params.sall_lla;
theta0 = measurement_params.theta0;
wE = measurement_params.wE;

% Initialize
xbar = xbar0;
batch_cnt = 0;
err = 0.1;
n = length(xbar0);
Xhat(:,:,1) = Xnom;

% Iterative algorithm
while err > tol
    k = batch_cnt + 1;
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
    [~,sstm] = ode45(@(t,x) odeSTM_J2_rv(t,x,mu,Rp,J2),t,state_stm,options);
    Xhat(:,:,k) = sstm(:,1:n);

    % Generate measurements for this initial conditions
    [G,~] = genMeasurements(t,stations_lla,theta0,wE,[],sstm(:,1:n));

    for i = 1:length(t)

        % Get this time step's state info
        Xbari(i,:) = sstm(i,1:n);
        Phii(:,:,i) = reshape(sstm(i,n+1:end),n,n);
    
        % Read next observation and expected observation
        M = Y{i};

        % Computes observation deviation
        if isempty(M)
            y(:,i,k) = NaN;
            continue % skips measurement update if no measurement available
        end
        
        % Select correct station measurement
        Mn = G{i};
        for j = 1:size(Mn,1)
            if Mn(j,1) ~= M(1)
                continue
            else
                Mnj = Mn(j,:);
                break
            end
        end

        % Prefits and update
        y(:,i,k) = M(:,2:3)' - Mnj(:,2:3)';
        Htilde = sc_range_ranger_Htilde(Xbari(i,:),gs_state(i,:));
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
    Phi_i0 = Phii(:,:,1);
    for i = 1:length(t)
        % if i > 1
        %     Phi_i0 = Phii(:,:,i) * Phi_i0;   % Φ(ti,t0) = Φ(ti,ti-1) * Φ(ti-1,t0)
        % end
        P(:,:,i,k) = Phii(:,:,i)*P0*Phii(:,:,i)';
        Htilde = sc_range_ranger_Htilde(Xbari(i,:),gs_state(i,:));
        Hi = Htilde*Phii(:,:,i);
        yhat(:,i,k) = y(:,i,k) - Hi*dxhat;
    end

    % Increase counter
    batch_cnt = batch_cnt + 1;
end

% Transpose outputs
yhat = pagetranspose(yhat);

end