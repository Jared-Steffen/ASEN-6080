function [cov_BPlane,BdotRhat,BdotThat] = BPlaneAnalysisV3(DCO_t,N,DCO_state,DCO_covariance,constants)

% Set RSOI Gate
RSOI_gate = N*925000;

% Getinitial conditions
t0 = DCO_t(end);
X0 = DCO_state(end,:);
P0 = DCO_covariance(:,:,end);
n = length(X0);

% Set ode45 options
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events', @(t,s) BPlaneEvents(t,s,RSOI_gate));

% Integrate to RSOI_gate
tspan = [t0 t0+1e8]; % Really large upper bound, won't be hit
[tint, Xint, te, Xe, ie] = ode45(@(t,x) odeSTM_SunEarth3BP_SRP3(t,x,constants), tspan, [X0(:);reshape(eye(n),[],1)], options);

% Extract state
Xf = Xe(end,1:n);
tf = te(end);
r = Xf(1:3);
v = Xf(4:6);


% Get eci2str DCM
muE = constants.muEarth;
e_vec = 1/muE.*((norm(v)^2-muE/norm(r)).*r(:) - (dot(r(:),v(:))).*v(:));
e = norm(e_vec);
Phat = e_vec./norm(e_vec);
h_vec = cross(r(:),v(:));
What = h_vec./norm(h_vec);
Qhat = cross(What,Phat);
a = -muE/(norm(v)^2-(2*muE)/norm(r));
p = norm(h_vec)^2/muE;
% e = 1-p/a;
b = norm(a)*sqrt(e^2-1);
Shat = v(:)/norm(v);
Nhat = [0 0 1]';
That = cross(Shat,Nhat)./norm(cross(Shat,Nhat));
Rhat = cross(Shat,That);
B = b.*cross(Shat,What);

eci2str = [Shat';That';Rhat'];

% Get LTOF
nu = acos(dot(r/norm(r),Phat));
f = acosh(1+norm(v)^2/muE*(a*(1-e^2))/(1+e*cos(nu)));
LTOF = muE/norm(v)^3*(sinh(f)-f);

% Reintegrate from DCO to Bplane
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
tspan = [t0 tf+LTOF];
[~, Xint] = ode45(@(t,x) odeSTM_SunEarth3BP_SRP3(t,x,constants), tspan, [X0(:);reshape(eye(n),[],1)], options);

% STM
flatSTM = Xint(end,n+1:end);
Phi = reshape(flatSTM(:),n,n,[]);
PBplane = Phi*P0*Phi';

% Outputs
cov_BPlane = eci2str*PBplane(1:3,1:3)*eci2str';
BdotRhat = dot(B,Rhat);
BdotThat = dot(B,That);

end

