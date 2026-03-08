function [Xhat,P,y,yhat] = SRIF(t,xbar0,Pbar0,Q,Y,V,Xnom,constants,stations)
%{
SMOOTHER STILL TO BE FULLY IMPLEMENTED
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Q: process noise covariance matrix
    >Y: simulated actual measurement data for true trajectory
    >V: square root of measurement covariance matrix (upper triangular)
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
    >Xhat: full predicted state (dxhat + Xnom)/(used in smoother) 
    >dxhat: state estimate correction (used in smoother)
    >P: state history of estimated state covariance
    >y: pre-fit measurement residuals
    >yhat: post-fit residuals
    >Phii: STM from (t_k-1, tk) for each time step (used in smoother)
%}

% Set ode tolerance
options = odeset('RelTol',1e-11,'AbsTol',1e-11);

% Extract info
station_ids = stations.station_ids;

% Initialize
n = length(xbar0);
P = zeros(n,n,length(t));
y = NaN(2, length(t));
yhat = NaN(2, length(t));
X0 = Xnom(1,1:n);  
Lambda0 = inv(Pbar0);
R0 = chol(Lambda0);
Rim1 = R0;
xhatim1 = xbar0;
bim1 = Rim1*xhatim1;
Ru = chol(inv(Q));
bbaru = zeros(3,1);

% Integrate current initial condition
[~,Xint] = ode45(@(t,x) odeSTM_J2(t,x,constants),t,[X0';reshape(eye(n),[],1)],options);
Xbar = Xint(:,1:n);

% Iterative algorithm
for i = 1:length(t)

    % Get this time step's state info
    Xbari = Xbar(i,:); % Used for Htilde matrix
    Phi(:,:,i) = reshape(Xint(i,n+1:end),n,n);

    % Read next observation and determine station to pass to measurements
    M = Y(i,:);
    current_id = M(2);
    idx = find(current_id == station_ids);
    current_station = [current_id, stations.Rs(idx,:)];

    % Time update
    if i > 1
        Phii(:,:,i) = Phi(:,:,i)/Phi(:,:,i-1);
        deltat = t(i) - t(i-1);
    else
        Phii(:,:,i) = Phi(:,:,i);
        deltat = 0;
    end

    xbari = Phii(:,:,i)*xhatim1;

    % If adding process noise, then implement SNC
    if ~isempty(Q)  
        
        % Skip injection of SNC if delta t is big
        if deltat > 10
            Rbari = Rim1/Phii(:,:,i);
            bbari = Rbari*xbari;
    
            % Upper triangularize
            Rbari_size = size(Rbari,1);
            R_bbari = qr([Rbari bbari(:)]);
            Rbari = R_bbari(1:Rbari_size,1:Rbari_size);
            bbari = R_bbari(1:Rbari_size,end);
        else
            % SNC implementation
            Gamma = deltat .* [deltat/2.*eye(3);
                                      eye(3)];
            Rtilde = Rim1/Phii(:,:,i);
            
            % Householder transformation
            [Rbaru(:,:,i),Rbarux(:,:,i),butilde(:,i),Rbari,bbari] = householderTimeUpdate(Ru,bbaru,Rtilde,Gamma,bim1);
            ubar = Rbaru(:,:,i)\(butilde( :,i) - Rbarux(:,:,i)*xbari);
            bbaru = Ru*ubar;

        end
    else
        Rbari = Rim1/Phii(:,:,i);
        bbari = Rbari*xbari;

        % Upper triangularize
        Rbari_size = size(Rbari,1);
        R_bbari = qr([Rbari bbari(:)]);
        Rbari = R_bbari(1:Rbari_size,1:Rbari_size);
        bbari = R_bbari(1:Rbari_size,end);
        
    end

    % Extract station id and generate measurement
    [G,gs_state] = genSingleMeasurement(t(i),current_station,constants,Xbari(1:6));
    y(:,i) = M(:,3:4)' - G(:,2:3)';
    Htilde = linearizedH(Xbari(1:6),gs_state);

    % Whiten measurement model
    [ytilde,Htilde2] = whitenMeasModel(y(:,i),V,Htilde);

    % Householder transformation (measurement update)
    [R,b] = householderMeasUpdate(Rbari,bbari,Htilde2,ytilde);
    dxhat(:,i) = R\b;
    P(:,:,i) = inv(R)/R';
    if deltat > 10
        uhat(:,:,i) = zeros(3,1);
    else
        uhat(:,:,i) = Rbaru(:,:,i)\(butilde(:,i) - Rbarux(:,:,i)*dxhat(:,i));
    end

    % Post-fit residuals
    yhat(:,i) = y(:,i) - Htilde*dxhat(:,i);

    % Move iteration forward
    xhatim1 = dxhat(:,i);
    Rim1 = R;
    bim1 = b;
    
end

% Add state deviations to full state
Xhat = Xbar + dxhat(:,:)';

% Trasponse residuals
y = y';
yhat = yhat';

end

