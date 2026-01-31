function [Xhat,P,y,yhat] = orbitLKF(t,xbar0,Pbar0,Y,constants,measurement_params,R,Xnom,gs_state,num_iterations)
%{
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Y: simulated actual measurement data for true trajectory
    >constants: struct that contains:
        -mu: gravitational parameter for central body
        -Rp: radius of central body
        -J2: J2 coefficient of central body
    >measurement_params: struct containing measurement station locations in
                         lla, initial angle of central body spin, angular
                         velocity of central body, elevation mask in
                         degrees
    >R: measurement noise covariance matrix
    >Xnom: nominal trajectory state vector
    >gs_state: state of the ground stations in eci frame at each time step
    >num_iterations: number of times to iterate LKF algorithm

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
n = length(xbar0);
dxhat = zeros(n, length(t));
Phi = zeros(n, n, length(t));
P = zeros(n, n, length(t), num_iterations);
y = NaN(2, length(t), num_iterations);
yhat = NaN(2, length(t), num_iterations);
Xhat = zeros(length(t), n, num_iterations);
X0 = Xnom(1,:);

for j = 1:num_iterations
    
    % Reset for this iteration
    Pim1 = Pbar0;
    xhatim1 = xbar0;
    
    % Integrate current initial condition
    state_stm = [X0';reshape(eye(n),[],1)];
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
    [~,sstm] = ode45(@(t,x) odeSTM_J2_rv(t,x,mu,Rp,J2),t,state_stm,options);
    Xref = sstm(:,1:n);

    % Generate measurements for this initial conditions
    [G,~] = genMeasurements(t,stations_lla,theta0,wE,[],Xref);

    % Iterative algorithm
    for i = 1:length(t)
    
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
            P(:,:,i,j) = Pbari;
            y(:,i,j) = NaN;
            yhat(:,i,j) = NaN;
    
            % Move iteration forward
            xhatim1 = dxhat(:,i);
            Pim1 = P(:,:,i,j);
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
        y(:,i,j) = yi - Htilde*xbari;
        dxhat(:,i) = xbari + Ki*y(:,i,j);
        P(:,:,i,j) = (eye(n) - Ki*Htilde)*Pbari*(eye(n) - Ki*Htilde)' + Ki*R*Ki';
        yhat(:,i,j) = yi - Htilde*dxhat(:,i);
    
        % Move iteration forward
        xhatim1 = dxhat(:,i);
        Pim1 = P(:,:,i,j);
        
    end

    % Add deviations to nominal trajectory
    Xhat(:,:,j) = Xref + dxhat';

    % Backwards solve for new xhat0
    x_hat0 = Phi(:,:,end)\dxhat(:,end);
    X0 = X0 + x_hat0';
    xbar0 = xbar0 - x_hat0;
end


% Trasponse residuals
y = pagetranspose(y);
yhat = pagetranspose(yhat);

end

