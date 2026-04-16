function [Xhatm,dxhatm,Pm] = SRIFSmoother(t,dxhat,Xnom,P,Pbar,Phi,Gamma,btilde,Rbarux,Rbaru)
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
n = length(dxhat(:,end));
dxhatm = zeros(size(dxhat));
Pm = zeros(size(P));
dxhatm(:,end) = dxhat(:,end);
Pm(:,:,end) = P(:,:,end);

% Smoother algorithm
for i = length(t)-1:-1:1
    if btilde(:,i+1) == zeros(3,1)
        dxhatm(:,i) = Phi(:,:,i+1)\dxhatm(:,i+1);
        Sk = Phi(:,:,i+1)\eye(n);
        Pm(:,:,i) = P(:,:,i) + Sk*(Pm(:,:,i+1)-Pbar(:,:,i+1))*Sk';
    else
        uhatm = Rbaru(:,:,i+1)\(btilde(:,i+1)-Rbarux(:,:,i+1)*dxhatm(:,i+1));
        dxhatm(:,i) = Phi(:,:,i+1)\(dxhatm(:,i+1)-Gamma(:,:,i+1)*uhatm);
        Sk = Phi(:,:,i+1)\(eye(n)+Gamma(:,:,i+1)/Rbaru(:,:,i+1)*Rbarux(:,:,i+1));
        Pm(:,:,i) = P(:,:,i) + Sk*(Pm(:,:,i+1)-Pbar(:,:,i+1))*Sk';
    end
end

% Get full state
Xhatm = Xnom + dxhatm';


end