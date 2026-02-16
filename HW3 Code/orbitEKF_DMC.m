function [Xhat,P,y,yhat] = orbitEKF_DMC(t,xbar0,Pbar0,Q,Y,R,Xnom,constants,stations,LKFinit)
%{
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Q: process noise covariance components
    >Y: simulated actual measurement data for true trajectory
    >R: measurement noise covariance matrix
    >constants: struct that contains:
        -mu: Earth's gravitational parameter mu
        -J2: Earth's J2 coefficient
        -RE: Earth's radius
        -tau: time constant for DMC
    >stations: struct containing information on stations:
        -Rs: positions of each station in the ecef frame
        -station_ids: stations ids correspoding to stations in Rs
        -theta0: initial spin angle of Earth in simulation in radians
        -el_mask: elevation mask for GS to S/C visibility in radians
    >LKFinit: number of measurements to process with the LKF to initialize
Outputs:
    >Xhat: full predicted state (dxhat + Xnom)
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
beta = 1/constants.tau;
varQ = diag(Q);
Phi = zeros(n-3,n-3,length(t));
Pim1 = Pbar0;
xhatim1 = xbar0;
X0 = Xnom(1,:);
whatim1 = xbar0(7:end);

if LKFinit > 0
    % Determine time step to stop LKF initialization at
    tLKF = t(1:LKFinit);
    
    % Integrate current initial condition
    Xref = [X0';reshape(eye(n-3),[],1)];
    [~,Xint] = ode45(@(t,x) odeSTM_J2_DMC(t,x,constants),tLKF,Xref,options);
    Xbar = Xint(:,1:n);
    
    
    % Iterative algorithm
    for i = 1:length(tLKF)
    
        % Get this time step's state info
        Xbari = Xbar(i,:);
        Phi(:,:,i) = reshape(Xint(i,n+1:end),n-3,n-3);
    
        % Read next observation and determine station to pass to measurements
        M = Y(i,:);
        current_id = M(2);
        idx = find(current_id == station_ids);
        current_station = [current_id, stations.Rs(idx,:)];
    
        % Time update
        if i > 1
            Phii = Phi(:,:,i)/Phi(:,:,i-1);
            deltat = t(i) - t(i-1);
        else
            Phii = Phi(:,:,i);
            deltat = 0;
        end
        
        % DMC components of new STM
        wt = whatim1.*exp(-beta*deltat);
        vwt = whatim1./beta.*(1-exp(-beta*deltat));
        rwt = whatim1.*deltat/beta + whatim1./beta^2.*(exp(-beta*deltat)-1);
        phirw = diag(rwt./whatim1);
        phivw = diag(vwt./whatim1);
        phiww = diag(wt./whatim1);    
        Phii_aug(:,:,i) = [Phii, [phirw;phivw]; zeros(3,6), phiww];

        % Process noise covariance matrix    
        Qwrr = diag(varQ./beta^2.*(deltat^3/3-deltat^2/beta+deltat/beta^2*(1-2*exp(-beta*deltat))+1/(2*beta^3)*(1-exp(-2*beta*deltat))));
        Qwrv = diag(varQ./beta^2.*(0.5*deltat^2-deltat/beta*(1-exp(-beta*deltat))+1/beta^2*(0.5-exp(-beta*deltat)+0.5*exp(-2*beta*deltat))));
        Qwvv = diag(varQ./beta^2.*(deltat-1/beta*(1.5+0.5*exp(-2*beta*deltat)-2*exp(-beta*deltat))));
        Qwrw = diag(varQ./beta^2.*(1/(2*beta)*(1-exp(-2*beta*deltat))-deltat*exp(-beta*deltat)));
        Qwvw = diag(varQ./beta^2.*(0.5*(1+exp(-2*beta*deltat))-exp(-beta*deltat)));
        Qwww = diag(varQ./(2*beta).*(1-exp(-2*beta*deltat)));
        Qw = [Qwrr Qwrv Qwrw;
              Qwrv Qwvv Qwvw;
              Qwrw Qwvw Qwww];

        % Propagate
        xbari = Phii_aug(:,:,i)*xhatim1;
        Pbari = Phii_aug(:,:,i)*Pim1*Phii_aug(:,:,i)' + Qw;

        % Extract station id and generate measurement
        [G,gs_state] = genSingleMeasurement(t(i),current_station,constants,Xbari(1:6));

        y(:,i) = M(:,3:4)' - G(:,2:3)';
        Htilde = linearizedH_DMC(Xbari(1:6),gs_state);
        Ki = Pbari*Htilde'/(Htilde*Pbari*Htilde' + R);
    
        % Measurement correction
        dxhat(:,i) = xbari + Ki*(y(:,i)- Htilde*xbari);
        P(:,:,i) = (eye(n) - Ki*Htilde)*Pbari*(eye(n) - Ki*Htilde)' + Ki*R*Ki';
        yhat(:,i) = y(:,i) - Htilde*dxhat(:,i);
    
        % Move iteration forward
        whatim1 = dxhat(7:end,i);
        xhatim1 = dxhat(:,i);
        Pim1 = P(:,:,i);

    end
    
    % Add deviations to nominal trajectory for LKF section and initialize EKF
    Xhat = Xnom(1:length(tLKF),:)' + dxhat;
end



% Iterative algorithm for EKF
for i = LKFinit+1:length(t)

    if i == 1
        Xbar(i,:) = X0;
        Phii  = eye(n-3);
        deltat = 0;
        whatim1 = X0(7:end);
    else
        % Integrate reference trajectory and stm to next time step
        current_tstep = [t(i-1) t(i)];
        Xref = [Xhat(:,end);reshape(eye(n-3),[],1)];
    
        [~,Xint] = ode45(@(t,x) odeSTM_J2_DMC(t,x,constants),current_tstep,Xref,options);
        Xint_end = Xint(end,:);
        Xbar(i,:) = Xint_end(1:n);
        Phii = reshape(Xint_end(n+1:end),n-3,n-3);
        deltat = t(i) - t(i-1);
        whatim1 = Xhat(7:end,end);
    end

    % Read next observation and determine station to pass to measurements
    M = Y(i,:);
    current_id = M(2);
    idx = find(current_id == station_ids);
    current_station = [current_id, stations.Rs(idx,:)];

    % DMC components of new STM
    wt = whatim1.*exp(-beta*deltat);
    vwt = whatim1./beta.*(1-exp(-beta*deltat));
    rwt = whatim1.*deltat/beta + whatim1./beta^2.*(exp(-beta*deltat)-1);
    phirw = diag(rwt./whatim1);
    phivw = diag(vwt./whatim1);
    phiww = diag(wt./whatim1);    
    Phii_aug = [Phii, [phirw;phivw]; zeros(3,6), phiww];

    % Process noise covariance matrix    
    Qwrr = diag(varQ./beta^2.*(deltat^3/3-deltat^2/beta+deltat/beta^2*(1-2*exp(-beta*deltat))+1/(2*beta^3)*(1-exp(-2*beta*deltat))));
    Qwrv = diag(varQ./beta^2.*(0.5*deltat^2-deltat/beta*(1-exp(-beta*deltat))+1/beta^2*(0.5-exp(-beta*deltat)+0.5*exp(-2*beta*deltat))));
    Qwvv = diag(varQ./beta^2.*(deltat-1/beta*(1.5+0.5*exp(-2*beta*deltat)-2*exp(-beta*deltat))));
    Qwrw = diag(varQ./beta^2.*(1/(2*beta)*(1-exp(-2*beta*deltat))-deltat*exp(-beta*deltat)));
    Qwvw = diag(varQ./beta^2.*(0.5*(1+exp(-2*beta*deltat))-exp(-beta*deltat)));
    Qwww = diag(varQ./(2*beta).*(1-exp(-2*beta*deltat)));
    Qw = [Qwrr Qwrv Qwrw;
          Qwrv Qwvv Qwvw;
          Qwrw Qwvw Qwww];

    % Propagate
    Pbari = Phii_aug*Pim1*Phii_aug' + Qw;

    % Extract station id and generate measurement
    [G,gs_state] = genSingleMeasurement(t(i),current_station,constants,Xbar(i,1:6));

    y(:,i) = M(:,3:4)' - G(:,2:3)';
    Htilde = linearizedH_DMC(Xbar(i,1:6),gs_state);
    Ki = Pbari*Htilde'/(Htilde*Pbari*Htilde' + R);

    % Measurement correction
    dXhat = Ki*y(:,i);
    Xhat(:,i) = Xbar(i,:)' + dXhat;
    yhat(:,i) = y(:,i) - Htilde*dXhat;
    P(:,:,i) = (eye(n) - Ki*Htilde)*Pbari*(eye(n) - Ki*Htilde)' + Ki*R*Ki';

    % Move iteration forward
    whatim1 = Xhat(7:9,i);
    Pim1 = P(:,:,i);

end

% Transpose vectors
Xhat = Xhat';
yhat = yhat';
y = y';

end