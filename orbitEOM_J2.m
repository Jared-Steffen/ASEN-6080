function [var_dot] = orbitEOM_J2(~,var,mu,RE,J2)
    % Goal: Output ODEs for ode45

    % Extract state variables
    x = var(1);
    y = var(2);
    z = var(3);
    u = var(4);
    v = var(5);
    w = var(6);

    % Calculate radius
    r = sqrt(x^2 + y^2 + z^2);

    % J2 Acceleration Coefficient
    uJ2_coeff = (1.5*mu*RE^2*J2)/r^5;

    % Assign t.r.o.c variables
    x_dot = u;
    y_dot = v;
    z_dot = w;
    u_dot = (-mu*x)/r^3 + uJ2_coeff*(5*z^2/r^2-1)*x;
    v_dot = (-mu*y)/r^3 + uJ2_coeff*(5*z^2/r^2-1)*y;
    w_dot = (-mu*z)/r^3 + uJ2_coeff*(5*z^2/r^2-3)*z;

    % Final state derivative
    var_dot = [x_dot;y_dot;z_dot;u_dot;v_dot;w_dot;0];
end