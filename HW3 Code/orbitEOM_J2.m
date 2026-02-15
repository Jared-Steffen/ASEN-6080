function [var_dot] = orbitEOM_J2(~,var,constants)
    % Goal: Output ODEs for ode45

    % Extract constants
    mu = constants.mu;
    J2 = constants.J2;
    RE = constants.RE;

    % Extract state variables
    x = var(1);
    y = var(2);
    z = var(3);
    vx = var(4);
    vy = var(5);
    vz = var(6);

    % Calculate radius
    r = sqrt(x^2 + y^2 + z^2);

    % J2 Acceleration Coefficient
    uJ2_coeff = (1.5*mu*RE^2*J2)/r^5;

    % Assign t.r.o.c variables
    x_dot = vx;
    y_dot = vy;
    z_dot = vz;
    vx_dot = (-mu*x)/r^3 + uJ2_coeff*(5*z^2/r^2-1)*x;
    vy_dot = (-mu*y)/r^3 + uJ2_coeff*(5*z^2/r^2-1)*y;
    vz_dot = (-mu*z)/r^3 + uJ2_coeff*(5*z^2/r^2-3)*z;

    % Final state derivative
    var_dot = [x_dot;y_dot;z_dot;vx_dot;vy_dot;vz_dot];
end