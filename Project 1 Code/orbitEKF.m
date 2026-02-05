function [Xhat,P,y,yhat] = orbitEKF(t,xbar0,Pbar0,Y,R,Xnom,gs_state,constants,measurement_params)
%{
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Y: simulated actual measurement data for true trajectory
    >R: measurement noise covariance matrix
    >gs_state: state of the ground stations in eci frame at each time step
    >constants: struct that contains:
        -mu: gravitational parameter for central body
        -Rp: radius of central body
        -J2: J2 coefficient of central body
    >measurement_params: struct containing measurement station locations in
                         lla, initial angle of central body spin, angular
                         velocity of central body, elevation mask in
                         degrees
Outputs:
    >Xhat: full predicted state (dxhat + Xnom)
    >P: state history of estimated state covariance
    >y: pre-fit measurement residuals
    >yhat: post-fit residuals
%}

% Extract measurement params
stations_lla = measurement_params.sall_lla;
theta0 = measurement_params.theta0;
wE = measurement_params.wE;

% Extract constatns
mu = constants.mu;
J2 = constants.J2;
Rp = constants.Rp;

% Initialize
Pim1 = Pbar0;
xhatim1 = xbar0;
n = length(xbar0);
X0 = Xnom(1,:);

% Determine time step to stop LKF initialization at
obs_counter = 0;
for i = 1:length(t)
    M = Y{i};
        if isempty(M)
            continue
        end
        obs_counter = obs_counter + 1;
        if obs_counter == 100
            tLKF = t(1:i);
            tEKF = t(i:end);
            break
        end
end     

% Iterative algorithm for LKF initialization
state_stm = [X0';reshape(eye(n),[],1)];
options = odeset('RelTol',1e-11,'AbsTol',1e-11);
[~,sstm] = ode45(@(t,x) odeSTM_J2_rv(t,x,mu,Rp,J2),t,state_stm,options);
Xref = sstm(:,1:n);

% Generate measurements for this initial conditions
[G,~] = genMeasurements(t,stations_lla,theta0,wE,[],Xref);

    % Iterative algorithm
    for i = 1:length(tLKF)
    
        % Read next observation and expected observation
        M = Y{i};
        Mn = G{i};
    
        % Time update
        Phi(:,:,i) = reshape(sstm(i,n+1:end),n,n);
        if i > 1
            Phii = Phi(:,:,i)/Phi(:,:,i-1);
        else
            Phii = Phi(:,:,i);
        end
        xbari = Phii*xhatim1;
        Pbari = Phii*Pim1*Phii';
    
        % Computes observation deviation and Kalman Gain
        if isempty(M)
            dxhat(:,i) = xbari;
            P(:,:,i) = Pbari;
            y(:,i) = NaN;
            yhat(:,i) = NaN;
    
            % Move iteration forward
            xhatim1 = dxhat(:,i);
            Pim1 = P(:,:,i);
            continue % skips measurement update if no measurement available
        end
        % Select correct station measurement
        for k = 1:size(Mn,1)
            if Mn(k,1) ~= M(1)
                continue
            else
                Mnj = Mn(k,:);
                break
            end
        end

        yi = M(:,2:3)' - Mnj(:,2:3)';
        Htilde = sc_range_ranger_Htilde(Xref(i,:),gs_state(i,:));
        Ki = Pbari*Htilde'/(Htilde*Pbari*Htilde' + R);
    
        % Measurement correction
        y(:,i) = yi - Htilde*xbari;
        dxhat(:,i) = xbari + Ki*y(:,i);
        P(:,:,i) = (eye(n) - Ki*Htilde)*Pbari*(eye(n) - Ki*Htilde)' + Ki*R*Ki';
        yhat(:,i) = yi - Htilde*dxhat(:,i);
    
        % Move iteration forward
        xhatim1 = dxhat(:,i);
        Pim1 = P(:,:,i);
        
    end

% Add deviations to nominal trajectory for LKF section and initialize EKF
Xhat = Xnom(1:length(tLKF),:)' + dxhat;
Xhatim1 = Xhat(:,end)'; % different than xhatim1


% Iterative algorithm for EKF
for i = length(tLKF)+1:length(tLKF) + length(tEKF) - 1

    % Read next observation and expected observation
    M = Y{i};

    % Integrate reference trajectory and stm to next time step
    current_tstep = [t(i-1) t(i)];
    state_stm = [Xhatim1(:);reshape(eye(n),[],1)];
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
    [~,sstm] = ode45(@(t,x) odeSTM_J2_rv(t,x,mu,Rp,J2),current_tstep,state_stm,options);
    last_sstm = sstm(end,:);
    Xbari = last_sstm(1:n);
    Phii = reshape(last_sstm(n+1:end),n,n);

    % Time update
    Pbari = Phii*Pim1*Phii';

    % Computes observation deviation and Kalman Gain
    if isempty(M)
        Xhat(:,i) = Xbari;
        P(:,:,i) = Pbari;
        y(:,i) = NaN;
        yhat(:,i) = NaN;

        % Move iteration forward
        Xhatim1 = Xbari;
        Pim1 = Pbari;
        continue % skips measurement update if no measurement available
    end
    [G,~] = genMeasurements(t(i),stations_lla,theta0,wE,[],Xbari);
    Mn = G{1};
    % Select correct station measurement
    for k = 1:size(Mn,1)
        if Mn(k,1) ~= M(1)
            continue
        else
            Mnj = Mn(k,:);
            break
        end
    end
    y(:,i) = M(:,2:3)' - Mnj(:,2:3)';
    Htilde = sc_range_ranger_Htilde(Xbari,gs_state(i,:));
    Ki = Pbari*Htilde'/(Htilde*Pbari*Htilde' + R);

    % Measurement correction
    dXhat = Ki*y(:,i);
    Xhat(:,i) = Xbari(:) + dXhat;
    yhat(:,i) = y(:,i) - Htilde*dXhat;
    P(:,:,i) = (eye(n) - Ki*Htilde)*Pbari;

    % Move iteration forward
    Xhatim1 = Xhat(:,i)';
    Pim1 = P(:,:,i);

end

% Transpose vectors
Xhat = Xhat';
yhat = yhat';
y = y';

end
