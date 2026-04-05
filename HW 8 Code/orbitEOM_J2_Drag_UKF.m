function [var_dot] = orbitEOM_J2_Drag_UKF(~,var,constants)
    % Goal: Output ODEs for ode45

    % Extract constants
    RE = constants.RE;
    wE = constants.wE;
    rho0 = constants.rho0;
    r0 = constants.r0;
    H = constants.H;
    mu = constants.mu;
    J2 = constants.J2;
    A = constants.A;
    m = constants.m;
    CD = constants.CD;
    
    % Extract state variables
    x = var(1:6:end);
    y = var(2:6:end);
    z = var(3:6:end);
    vx = var(4:6:end);
    vy = var(5:6:end);
    vz = var(6:6:end);

    % Calculate radius
    r = sqrt(x.^2 + y.^2 + z.^2);
    vrel = sqrt((vx+wE.*y).^2+(vy-wE.*x).^2+vz.^2);

    % Calculate air density rho
    rho_atm = rho0*exp(-(r-r0)/H);

    % J2 Acceleration Coefficient
    J2_term = (1.5*mu*RE^2*J2)./r.^5;
    D_term = (rho_atm*CD*A.*vrel)./(2*m);

    % Assign t.r.o.c variables
    x_dot = vx;
    y_dot = vy;
    z_dot = vz;
    vx_dot = (-mu.*x)./r.^3 + J2_term.*(5.*z.^2./r.^2-1).*x - D_term.*(vx+wE.*y);
    vy_dot = (-mu.*y)./r.^3 + J2_term.*(5.*z.^2./r.^2-1).*y - D_term.*(vy-wE.*x);
    vz_dot = (-mu.*z)./r.^3 + J2_term.*(5.*z.^2./r.^2-3).*z - D_term.*vz;

    % Final state derivative
    var_dot = reshape([x_dot(:)';y_dot(:)';z_dot(:)';vx_dot(:)';vy_dot(:)';vz_dot(:)'],[],1);
end