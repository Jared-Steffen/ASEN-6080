function plot_covariance_ellipsoids(P,filter_type)
% plot_covariance_ellipsoids
%{
Inputs:
    >Xhat: Nx6 estimated state history [km, km/s]
    >P: 6x6xN covariance history
    >filter_type: string label (e.g. "Batch", "CKF")
Outputs:
    > 3σ POSITION and VELOCITY axis-aligned covariance ellipsoids at final time
%}

%% Final state and covariance
PF = P(:,:,end);     


Pr = PF(1:3,1:3);   % symmetrize
Pv = PF(4:6,4:6);

nsig = 3;

%% Axis-aligned 3σ radii
sig_r = sqrt(diag(Pr));
sig_v = sqrt(diag(Pv));

rx = nsig*sig_r(1); ry = nsig*sig_r(2); rz = nsig*sig_r(3);
rvx = nsig*sig_v(1); rvy = nsig*sig_v(2); rvz = nsig*sig_v(3);

%% ===== Position ellipsoid =====
figure(); hold on; grid on; axis equal; view(3)

ellipsoid(0, 0, 0, rx, ry, rz)
% xlim(rF(1) + [-rx rx])
% ylim(rF(2) + [-ry ry])
% zlim(rF(3) + [-rz rz])

xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]')
title(filter_type + " Final-Time Position Covariance Ellipsoid (±3σ)")

%% ===== Velocity ellipsoid =====
figure(); hold on; grid on; axis equal; view(3)

ellipsoid(0, 0, 0, rvx, rvy, rvz)

% xlim(vF(1) + [-rvx rvx])
% ylim(vF(2) + [-rvy rvy])
% zlim(vF(3) + [-rvz rvz])

xlabel('v_x [km/s]'); ylabel('v_y [km/s]'); zlabel('v_z [km/s]')
title(filter_type + " Final-Time Velocity Covariance Ellipsoid (±3σ)")

end