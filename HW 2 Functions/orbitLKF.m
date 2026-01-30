function [dxhat,Xhat,P,y,yhat] = orbitLKF(t,xbar0,Pbar0,Phi,Y,Ynom,R,Xnom,gs_state)
%{
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Phi: state transition matrix (computed offline) for nominal trajectory
    >Y: simulated actual measurement data for true trajectory
    >Ynom: expected measurements about the nominal trajectory
    >R: measurement noise covariance matrix
    >Xnom: nominal trajectory state vector
    >gs_state: state of the ground stations in eci frame at each time step
Outputs:
    >dxhat: time history of estimated state correction
    >Xhat: full predicted state (dxhat + Xnom)
    >P: state history of estimated state covariance
    >y: pre-fit measurement residuals
    >yhat: post-fit residuals
%}

% Initialize
Pim1 = Pbar0;
xhatim1 = xbar0;
n = length(xbar0);
dxhat = zeros(n,length(t));
P     = zeros(n,n,length(t));
y     = zeros(2,length(t));
yhat  = zeros(2,length(t));

% Iterative algorithm
for i = 1:length(t)

    % Read next observation and expected observation
    M = Y{i};
    Mnom = Ynom{i};

    % Time update
    if i > 1
        Phii = Phi(:,:,i)/Phi(:,:,i-1);
    else
        Phii = eye(n);
    end
    xbari = Phii * xhatim1;
    Pbari = Phii*Pim1*Phii';

    % Computes observation deviation and Kalman Gain
    if isempty(M) || isempty(Mnom)
        dxhat(:,i) = xbari;
        P(:,:,i) = Pbari;
        y(:,i) = NaN;
        yhat(:,i) = NaN;

        % Move iteration forward
        xhatim1 = dxhat(:,i);
        Pim1 = P(:,:,i);
        continue % skips measurement update if no measurement available
    end
    yi = M(:,2:3)' - Mnom(:,2:3)';
    Htilde = sc_range_ranger_Htilde(Xnom(i,:),gs_state(i,:));
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

% Add deviations to nominal trajectory
Xhat = Xnom + dxhat';

% Trasponse residuals
y = y'; yhat = yhat';

end

