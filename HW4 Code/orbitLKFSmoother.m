function [Xhatl,Pl] = orbitLKFSmoother(t,dxhat,Xnom,P,Pbar,Phi)
%{
Inputs:
    >t: time vector
    >dxhat: state estimate correction from LKF
    >Xnom: nominal trajectory state vector
    >P: state history of estimated state covariance from LKF
    >Pbar: time update covariances from LKF
    >Phi: STM from (t_k-1, tk) for each time step from LKF

Outputs:
    >Xhatl: smoother full predicted state (dxhat + Xnom)
    >Pl: smoother state history of estimated state covariance
%}

% Initialize
dxhatl = zeros(size(dxhat));
Pl = zeros(size(P));
dxhatl(:,end) = dxhat(:,end);
Pl(:,:,end) = P(:,:,end);


% Smoother algorithm
for i = length(t)-1:-1:1
    Sk = P(:,:,i)*Phi(:,:,i+1)'/Pbar(:,:,i+1);
    dxhatl(:,i) = dxhat(:,i) + Sk*(dxhatl(:,i+1) - Phi(:,:,i+1)*dxhat(:,i));
    Pl(:,:,i) = P(:,:,i) + Sk*(Pl(:,:,i+1)-Pbar(:,:,i+1))*Sk';
end

% Get full state
Xhatl = Xnom + dxhatl';


end

