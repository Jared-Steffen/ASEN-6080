function [dxhat,Xhat,P,y] = orbitLKF(t,xbar0,Pbar0,Phi,Y,R,Xnom)
%{
Inputs:
    >t: time vector
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Phi: state transition matrix (computed offline) for nominal trajectory
    >Y: simulated measurement data for true trajectory
    >R: measurement noise covariance matrix
    >Xnom: nominal trajectory state vector
Outputs:
    >dxhat: time hsitory of estimated state correction
    >Xhat: full predicted state (dxhat + Xnom)
    >P: state history of estimated state covariance
    >y: measurement residuals
%}

% Initialize

% Time update

% Measurement Correction


end

