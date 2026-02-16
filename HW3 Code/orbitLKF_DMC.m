function [Xhat,P,y,yhat] = orbitLKF_DMC(t,xbar0,Pbar0,Q,Y,R,Xnom,constants,stations,num_iterations)
%{
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Q: process noise covariance components 
    >Y: simulated actual measurement data for true trajectory
    >R: measurement noise covariance matrix
    >Xnom: nominal trajectory state vector
    >constants: struct that contains:
        -mu: Earth's gravitational parameter mu
        -J2: Earth's J2 coefficient
        -RE: Earth's radius
        -wE: Earth's rotation rate
        -tau: time constant for DMC
    >stations: struct containing information on stations:
        -Rs: positions of each station in the ecef frame
        -station_ids: stations ids correspoding to stations in Rs
        -theta0: initial spin angle of Earth in simulation in radians
    >num_iterations: number of times to iterate LKF algorithm

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
P = zeros(n,n,length(t),num_iterations);
y = NaN(2, length(t),num_iterations);
yhat = NaN(2, length(t),num_iterations);
Xhat = zeros(length(t),n,num_iterations);
X0 = Xnom(1,1:n);

for j = 1:num_iterations
    
    % Reset for this iteration
    Pim1 = Pbar0;
    xhatim1 = xbar0;
    whatim1 = xbar0(7:end);

    % Integrate current initial condition
    [~,Xint] = ode45(@(t,x) odeSTM_J2_DMC(t,x,constants),t,[X0';reshape(eye(n-3),[],1)],options);
    Xbar = Xint(:,1:n);

    % Iterative algorithm
    for i = 1:length(t)

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

        y(:,i,j) = M(:,3:4)' - G(:,2:3)';
        Htilde = linearizedH_DMC(Xbari(1:6),gs_state);
        Ki = Pbari*Htilde'/(Htilde*Pbari*Htilde' + R);
    
        % Measurement correction
        dxhat(:,i) = xbari + Ki*(y(:,i,j)- Htilde*xbari);
        P(:,:,i,j) = (eye(n) - Ki*Htilde)*Pbari*(eye(n) - Ki*Htilde)' + Ki*R*Ki';
        yhat(:,i,j) = y(:,i,j) - Htilde*dxhat(:,i);
    
        % Move iteration forward
        whatim1 = dxhat(7:end,i);
        xhatim1 = dxhat(:,i);
        Pim1 = P(:,:,i,j);
        
    end

    % Add deviations to nominal trajectory
    Xhat(:,:,j) = Xbar + dxhat';

    % Backwards solve for new xhat0
    x_hat0 = Phii_aug(:,:,end)\dxhat(:,end);
    X0 = X0 + x_hat0';
    xbar0 = xbar0 - x_hat0;
end


% Trasponse residuals
y = pagetranspose(y);
yhat = pagetranspose(yhat);

end

