function [Xhat,dxhat,P,Pc,Pxx,S,y,yhat,Psii] = SRIFConsiderCov(t,xbar0,c,Pbar0,Pcc,Y,V,X0,constants,stations)
%{
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Q: process noise covariance matrix
    >Y: simulated actual measurement data for true trajectory
    >V: square root of measurement covariance matrix (upper triangular)
    >Xnom: nominal trajectory state vector
    >constants: struct that contains:
        -RE: Earth's radius
        -wE: Earth's rotation rate
    >stations: struct containing information on stations:
        -Rs: positions of each station in the ecef frame
        -station_ids: stations ids correspoding to stations in Rs
        -theta0: initial spin angle of Earth in simulation in radians

Outputs:
    >Xhat: full predicted state (dxhat + Xnom)/(used in smoother) 
    >dxhat: state estimate correction (used in smoother)
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
nc = length(c);
nt = length(t);
P = zeros(n,n,nt);
y = NaN(2, nt);
yhat = NaN(2, nt);
X0 = X0(1,1:n);
Sim1 = zeros(n,1);

% Initialize
Lambda0 = inv(Pbar0);
R0 = chol(Lambda0);
Rim1 = R0;
xhatim1 = xbar0;

% Integrate current initial condition
[~,Xint] = ode45(@(t,x) odeSTM_SunEarth3BP_SRPConsiderCov(t,x,constants),t,[X0';reshape(eye(n+nc),[],1)],options);
Xbar = Xint(:,1:n);

% Iterative algorithm
for i = 1:length(t)

    % Get this time step's state info
    Xbari = Xbar(i,:); % Used for Htilde matrix
    Psi(:,:,i) = reshape(Xint(i,n+1:end),n+nc,n+nc);

    % Read next observation and determine station to pass to measurements
    M = Y(i,:);
    current_id = M(2);
    idx = find(current_id == station_ids);
    current_station = [current_id, stations.Rs(idx,:)];

    % Time update
    if i > 1
        Psii(:,:,i) = Psi(:,:,i)/Psi(:,:,i-1);
        deltat = t(i) - t(i-1);
    else
        Psii(:,:,i) = Psi(:,:,i);
        deltat = 0;
    end
    
    % Extract necessary components
    Phii(:,:,i) = Psii(1:n,1:n,i);
    thetai(:,:,i) = Psii(1:n,n+1:end,i);

    % SRIF Time Update
    xbari = Phii(:,:,i)*xhatim1;
    Rbari = Rim1/Phii(:,:,i);
    bbari = Rbari*xbari;
    Rbari_size = size(Rbari,1);
    R_bbari = qr([Rbari bbari(:)]);
    Rbari = R_bbari(1:Rbari_size,1:Rbari_size);
    bbari = R_bbari(1:Rbari_size,end);
    Pbari = inv(Rbari)/Rbari';

    % Consider Time Update
    Sbari = Phii(:,:,i)*Sim1 + thetai(:,:,i);
    Pxxbari = Pbari + Sbari*Pcc*Sbari';
    Pxcbari = Sbari*Pcc;

    % Extract station id and generate measurement
    [G,gs_state] = genSingleMeasurement(t(i),current_station,constants,Xbari(1:6));
    y(:,i) = M(:,3:4)' - G(:,2:3)';
    Htilde = linearizedH(Xbari,gs_state);
    Hctilde = linearizedHc(t(i),Xbari,[current_id;gs_state],constants,station_ids);

    % Whiten measurement model
    [ytilde,Htilde2,Hctilde2] = whitenMeasModelConsiderCov(y(:,i),V,Htilde,Hctilde);

    % Householder transformation (measurement update)
    [R,b] = householderMeasUpdate(Rbari,bbari,Htilde2,ytilde);
    dxhat(:,i) = R\b;
    P(:,:,i) = inv(R)/R';

    % Consider measurement update
    Lambdabari = Rbari'*Rbari;
    q = Lambdabari*Sbari - Htilde2'*Hctilde2;
    Z = R'\q;
    S(:,:,i) = R\Z;
    Pxx(:,:,i) = P(:,:,i) + S(:,:,i)*Pcc*S(:,:,i)';
    Pxc(:,:,i) = S(:,:,i)*Pcc;

    % Post-fit residuals
    yhat(:,i) = y(:,i) - Htilde*dxhat(:,i);

    % Full consider cov matrix
    Pc(:,:,i) = [Pxx(:,:,i), Pxc(:,:,i);
                 Pxc(:,:,i)', Pcc];

    % Move iteration forward
    xhatim1 = dxhat(:,i);
    Rim1 = R;
    Sim1 = S(:,:,i);
    

end

% Add state deviations to full state
Xhat = Xbar + dxhat';

% Trasponse residuals
y = y';
yhat = yhat';

end
