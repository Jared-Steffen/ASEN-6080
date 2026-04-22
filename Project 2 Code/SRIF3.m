function [Xhat,Xhatm,dxhat,dxhatm,P,Pm,y,yhat,yhatl] = SRIF3(t,xbar0,Pbar0,Q,Y,V,X0,constants,stations,num_iterations)
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
    >num_iterations: number of times to iterate SRIF algorithm

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
warning('off','MATLAB:nearlySingularMatrix')

% Extract info
station_ids = stations.station_ids;

% Initialize
n = length(xbar0);
nt = length(t);

P = zeros(n,n,nt,num_iterations);
Pm = zeros(n,n,nt,num_iterations);
y = NaN(2, nt, num_iterations);
yhat = NaN(2, nt, num_iterations);
yhatl = NaN(2, nt, num_iterations);
dxhat = zeros(n,nt,num_iterations);
dxhatm = zeros(n,nt,num_iterations);
Xhat = zeros(nt,n,num_iterations);
Xhatm = zeros(nt,n,num_iterations);

X0 = X0(1,1:n);

for j = 1:num_iterations

    % Reset for this iteration
    L = chol(Pbar0, 'lower');
    R0 = L \ eye(size(L));
    Rim1 = R0;
    xhatim1 = xbar0;
    bim1 = Rim1*xhatim1;
    Ru = chol(inv(Q));
    bbaru = zeros(3,1);

    tm = constants.tm;

    T = find(t < tm, 1, 'last');
    
    x0_aug = [X0'; reshape(eye(n),[],1)];
    
    % Segment 1: propagate up to tm
    t1 = [t(1:T); tm];
    [~,Xint1] = ode45(@(t,x) odeSTM_SunEarth3BP_SRP3(t,x,constants), ...
                      t1, x0_aug, options);
    
    % Extract final state+STM at tm-
    x_tm_minus = Xint1(end,1:n).';
    Phi_tm_minus = reshape(Xint1(end,n+1:end),n,n);
    
    % Apply impulsive state jump at tm
    x_tm_plus = x_tm_minus;
    x_tm_plus(4:6) = x_tm_plus(4:6) + x_tm_plus(11:13);
    
    % Apply STM jump at tm
    J = eye(n);
    J(4:6,11:13) = eye(3);
    Phi_tm_plus = J * Phi_tm_minus;
    
    % Build new initial condition for segment 2
    x0_seg2 = [x_tm_plus; reshape(Phi_tm_plus,[],1)];
    
    % Segment 2: propagate from tm to end
    t2 = [tm; t(T+1:end)];
    [~,Xint2] = ode45(@(t,x) odeSTM_SunEarth3BP_SRP3(t,x,constants), ...
                      t2, x0_seg2, options);
    
    Xint1(end,:) = [];
    % Stitch results together without duplicating tm row
    Xint = [Xint1; Xint2(2:end,:)];
    
    % Nominal trajectory
    Xbar = Xint(:,1:n);


    % Iterative algorithm
    for i = 1:length(t)

        % Get this time step's state info
        Xbari = Xbar(i,:); % Used for Htilde matrix
        Phi(:,:,i,j) = reshape(Xint(i,n+1:end),n,n);

        % Read next observation and determine station to pass to measurements
        M = Y(i,:);
        current_id = M(2);
        idx = find(current_id == station_ids);
        if current_id == 65
            current_station = [current_id, Xbari(8:10)];
        else
            current_station = [current_id, stations.Rs(idx,:)];
        end

        % Time update
        if i > 1
            Phii(:,:,i,j) = Phi(:,:,i,j)/Phi(:,:,i-1,j);
            deltat = t(i) - t(i-1);
        else
            Phii(:,:,i,j) = Phi(:,:,i,j);
            deltat = 0;
        end

        xbari = Phii(:,:,i,j)*xhatim1;

        % If adding process noise, then implement SNC
        if ~isempty(Q)  
            
            % Skip injection of SNC if delta t is big
            if deltat > 60
                Rbari = Rim1/Phii(:,:,i,j);
                bbari = Rbari*xbari;

                % Upper triangularize
                Rbari_size = size(Rbari,1);
                R_bbari = qr([Rbari bbari(:)]);
                Rbari = R_bbari(1:Rbari_size,1:Rbari_size);
                bbari = R_bbari(1:Rbari_size,end);
            else
                % SNC implementation
                Gamma(:,:,i,j) = zeros(n,3);
                Gamma(1:6,:,i,j) = deltat .* [deltat/2.*eye(3);
                                              eye(3)];
                Rtilde = Rim1/Phii(:,:,i,j);
                
                % Householder transformation
                [Rbaru(:,:,i,j),Rbarux(:,:,i,j),butilde(:,i,j),Rbari,bbari] = ...
                    householderTimeUpdate(Ru,bbaru,Rtilde,Gamma(:,:,i,j),bim1);
                ubar = Rbaru(:,:,i,j)\(butilde(:,i,j) - Rbarux(:,:,i,j)*xbari);
                bbaru = Ru*ubar;
            end
        else
            Rbari = Rim1/Phii(:,:,i,j);
            bbari = Rbari*xbari;

            % Upper triangularize
            Rbari_size = size(Rbari,1);
            R_bbari = qr([Rbari bbari(:)]);
            Rbari = R_bbari(1:Rbari_size,1:Rbari_size);
            bbari = R_bbari(1:Rbari_size,end);
        end
        J = Rbari \ eye(size(Rbari));
        Pbari(:,:,i,j) = J*J';

        % Extract station id and generate measurement
        [G,gs_state] = genSingleMeasurement(t(i),current_station,constants,Xbari(1:6));
        y(:,i,j) = M(:,3:4)' - G(:,2:3)';
        Htilde(:,:,i,j) = linearizedH3(t(i),Xbari,[current_id;gs_state],constants,station_ids);

        % Whiten measurement model
        [ytilde,Htilde2] = whitenMeasModel(y(:,i,j),V,Htilde(:,:,i,j));

        % Householder transformation (measurement update)
        [R,b] = householderMeasUpdate(Rbari,bbari,Htilde2,ytilde);
        dxhat(:,i,j) = R\b;
        P(:,:,i,j) = inv(R)/R';
        if deltat > 60
            uhat(:,:,i,j) = zeros(3,1);
        else
            uhat(:,:,i,j) = Rbaru(:,:,i,j)\(butilde(:,i,j) - Rbarux(:,:,i,j)*dxhat(:,i,j));
        end

        % Post-fit residuals
        yhat(:,i,j) = y(:,i,j) - Htilde(:,:,i,j)*dxhat(:,i,j);

        % Move iteration forward
        xhatim1 = dxhat(:,i,j);
        Rim1 = R;
        bim1 = b;
        
    end

    [Xhatm(:,:,j),dxhatm(:,:,j),Pm(:,:,:,j)] = SRIFSmoother( ...
        t,dxhat(:,:,j),Xbar,P(:,:,:,j),Pbari(:,:,:,j),Phii(:,:,:,j), ...
        Gamma(:,:,:,j),butilde(:,:,j),Rbarux(:,:,:,j),Rbaru(:,:,:,j));

    % Get smoothed residuals
    for i = 1:length(t)
        yhatl(:,i,j) = y(:,i,j) - Htilde(:,:,i,j)*dxhatm(:,i,j);
    end

    % Add state deviations to full state
    Xhat(:,:,j) = Xbar + dxhat(:,:,j)';

    % Backwards solve for new xhat0
    x_hat0 = Phi(:,:,end,j)\dxhat(:,end,j);
    X0 = X0 + x_hat0';
    xbar0 = xbar0 - x_hat0;
end

% Trasponse residuals
y = pagetranspose(y);
yhat = pagetranspose(yhat);
yhatl = pagetranspose(yhatl);

end