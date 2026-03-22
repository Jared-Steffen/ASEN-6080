function [Xhat,Xchat,dxhat,dxchat,P,Pxx,Pc,Psi,y,yhat] = orbitLKF_consider(t,xbar0,c,Pbar0,Pcc,Y,R,Xnom,constants,stations)
%{
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >c: consider parameters
    >Pbar0: initialized state covariance matrix
    >Pxxbar0: initialized state covariance matrix including consider
    >Pcc: consider parameter covariance matrix
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

Outputs:
    >Xhat: full predicted state (dxhat + Xnom)
    >Xchat: full predicted state with consider (dxchat + Xnom)
    >dxhat: state estimate correction 
    >dxchat: state estimate correction  w/ consider parameter
    >P: history of state covariance
    >Pxx: history of state with consider covariance
    >Pc: full consider covariance matrix
    >Psi: full augmented STM
    >y: pre-fit measurement residuals
    >yhat: post-fit residuals
%}

% Set ode tolerance
options = odeset('RelTol',1e-11,'AbsTol',1e-11);

% Extract info
station_ids = stations.station_ids;

% Initialize
n = length(xbar0);
Psi = zeros(n+length(c),n+length(c),length(t));
P = zeros(n,n,length(t));
y = NaN(2, length(t));
yhat = NaN(2, length(t));
Xhat = zeros(length(t),n);
X0 = Xnom(1,1:n);
Pim1 = Pbar0;
Sim1 = zeros(n,1);
xhatim1 = xbar0;

% Integrate current initial condition
[~,Xint] = ode45(@(t,x) odeSTM_J2considerJ3(t,x,constants),t,[X0(:);reshape(eye(n+1),[],1)],options);
Xbar = Xint(:,1:n);

% Iterative algorithm
for i = 1:length(t)

    % Get this time step's state info
    Xbari = Xbar(i,:);
    Psi(:,:,i) = reshape(Xint(i,n+1:end),n+length(c),n+length(c));

    % Read next observation and determine station to pass to measurements
    M = Y(i,:);
    current_id = M(2);
    idx = find(current_id == station_ids);
    current_station = [current_id, stations.Rs(idx,:)];

    % Time update
    if i > 1
        Psii = Psi(:,:,i)/Psi(:,:,i-1);
        deltat = t(i) - t(i-1);
    else
        Psii = Psi(:,:,i);
        deltat = 0;
    end
    Phii = Psii(1:n,1:n);
    thetai = Psii(1:n,n+1);
    Pbari = Phii*Pim1*Phii';
    xbari = Phii*xhatim1;

    % Consider time update
    Sibar = Phii*Sim1 + thetai(:);
    xcbari = xbari + Sibar*c;
    Pxxbari = Pbari+ Sibar*Pcc*Sibar';
    Pxcbari = Sibar*Pcc;

    % Extract station id and generate measurement
    [G,gs_state] = genSingleMeasurement(t(i),current_station,constants,Xbari(1:6));

    y(:,i) = M(:,3:4)' - G(:,2:3)';
    Htilde = linearizedH(Xbari(1:6),gs_state);
    Hctilde = zeros(2,1); % J3 does not change measurement equations
    Ki = Pbari*Htilde'/(Htilde*Pbari*Htilde' + R);

    % Measurement correction
    dxhat(:,i) = xbari + Ki*(y(:,i)- Htilde*xbari);
    P(:,:,i) = (eye(n) - Ki*Htilde)*Pbari*(eye(n) - Ki*Htilde)' + Ki*R*Ki';
    yhat(:,i) = y(:,i) - Htilde*dxhat(:,i);

    % Consider measurement correction
    Si = (eye(n) - Ki*Htilde)*Sibar - Ki*Hctilde;
    dxchat(:,i) = dxhat(:,i) + Si*c;
    Pxx(:,:,i) = P(:,:,i) + Si*Pcc*Si';
    Pxc = Si*Pcc;

    % Form full consider covariance matrix
    Pc(:,:,i) = [Pxx(:,:,i) Pxc;
                 Pxc' Pcc];

    % Move iteration forward
    xhatim1 = dxhat(:,i);
    Pim1 = P(:,:,i);
    Sim1 = Si;
    
end

% Add deviations to nominal trajectory
Xhat(:,:) = Xbar + dxhat(:,:)';
Xchat(:,:) = Xbar + dxchat(:,:)';



% Trasponse residuals
y = y';
yhat = yhat';

end
