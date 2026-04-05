function [Xhat,P,y,yhat] = UKF_J2(t,X0,Pbar0,Q,alpha,Y,R,constants,stations)
%{
Inputs:
    >t: time vector
    >X0: initailized state
    >Pbar0: initialized state covariance matrix
    >Q: process noise covariance matrix
    >alpha: UKF tuning parameter
    >Y: simulated actual measurement data for true trajectory
    >R: measurement noise covariance matrix
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
    >Xhat: full predicted state
    >P: state history of estimated state covariance
    >y: pre-fit measurement residuals
    >yhat: post-fit residuals
%}

% Set ode tolerance
options = odeset('RelTol',1e-11,'AbsTol',1e-11);

% Extract info
station_ids = stations.station_ids;

% Initialize
n = length(X0);
Pim1 = Pbar0;
Xim1 = X0(:);

% Form constants and weights
kappa = 3-n;
beta = 2;
lambda = alpha^2*(n+kappa)-n;
gamma = sqrt(n+lambda);
Wm0 = lambda/(n+lambda);
Wc0 = lambda/(n+lambda)+(1-alpha^2+beta);
Wmc = 1/(2*(n+lambda));
Wm = [Wm0;Wmc*ones(2*n,1)];
Wc = [Wc0;Wmc*ones(2*n,1)];

% Iterative algorithm for UKF
for i = 1:length(t)

    if i == 1
        deltat = 0;
        S = chol(Pim1,'lower');
        Chii = [Xim1 Xim1+gamma.*S Xim1-gamma.*S];
    else
        % Integrate reference trajectory and stm to next time step
        current_tstep = [t(i-1) t(i)];
        S = chol(Pim1,'lower');
        Chiim1 = reshape([Xim1 Xim1+gamma.*S Xim1-gamma.*S],[],1);    
        [~,Chii_int] = ode45(@(t,x) orbitEOM_J2_UKF(t,x,constants),current_tstep,Chiim1,options);
        Chii = Chii_int(end,:);
        Chii = reshape(Chii,n,2*n+1);
        deltat = t(i) - t(i-1);
    end

    % Read next observation and determine station to pass to measurements
    M = Y(i,:);
    current_id = M(2);
    idx = find(current_id == station_ids);
    current_station = [current_id, stations.Rs(idx,:)];

    % Get new mean
    sum_m = zeros(size(X0(:)));
    for j = 1:size(Chii,2)
        sum_m = sum_m + Wm(j).*Chii(:,j); 
    end
    Xbari = sum_m;
    
    % Get new covariance
    sum_c = zeros(size(Pbar0));
    for j = 1:size(Chii,2)
        sum_c = sum_c + Wc(j).*(Chii(:,j)-Xbari)*(Chii(:,j)-Xbari)';
    end

    % SNC implementation
    if ~isempty(Q)
        Gamma = deltat .* [deltat/2.*eye(3);
                            eye(3)];
    
        % Skip SNC if measurement gap is big enough
        if deltat > 10
            Pbari = sum_c;
        else
            Pbari = sum_c + Gamma*Q*Gamma';
        end
    else
        Pbari = sum_c;
    end

    % New sigma points
    S = chol(Pbari,'lower');
    Chii = [Xbari Xbari+gamma.*S Xbari-gamma.*S];  

    % Push sigma points through nonlinear measurement model
    G = genMeasurementUKF(t(i),current_station,constants,Chii);

    % Get expected measurement
    sum_m = zeros(size(G,1),1);
    for j = 1:size(Chii,2)
        sum_m = sum_m + Wm(j).*G(:,j); 
    end
    ybar = sum_m;

    % Compute covariances
    sum_cy = zeros(size(R));
    sum_cxy = zeros(n,size(R,1));
    for j = 1:size(Chii,2)
        sum_cy = sum_cy + Wc(j).*(G(:,j)-ybar)*(G(:,j)-ybar)';
        sum_cxy = sum_cxy + Wc(j).*(Chii(:,j)-Xbari)*(G(:,j)-ybar)';
    end
    Pyy = sum_cy + R;
    Pxy = sum_cxy;

    % Measurement update
    Ki = Pxy*inv(Pyy);
    y(:,i) = M(:,3:4)'-ybar;
    Xhat(:,i) = Xbari + Ki*y(:,i);
    P(:,:,i) = Pbari - Ki*Pyy*Ki';

    % Calculate postfits
    Ghat = genMeasurementUKF(t(i),current_station,constants,Xhat(:,i));
    yhat(:,i) = M(:,3:4)' - Ghat;

    % Move iteration forward
    Pim1 = P(:,:,i);
    Xim1 = Xhat(:,i);

end

% Transpose vectors
Xhat = Xhat';
yhat = yhat';
y = y';


end