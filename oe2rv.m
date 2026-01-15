function [r,v] = oe2rv(mu,a,e,Omega,inc,w,nu)
    % Goal: Determine r and v at every time step of the simulation

    % Unit vectors
    x_hat = [1 0 0]';
    y_hat = [0 1 0]';
    z_hat = [0 0 1]';
    n_Omega_hat = cos(Omega)*x_hat + sin(Omega)*y_hat;
    n_Omega_hat_perp = -cos(inc)*sin(Omega)*x_hat + ...
        cos(inc)*cos(Omega)*y_hat + sin(inc)*z_hat;
    e_hat = cos(w)*n_Omega_hat + sin(w)*n_Omega_hat_perp;
    e_hat_perp = -sin(w)*n_Omega_hat + cos(w)*n_Omega_hat_perp;

    % Initialize r and v vectors
    r = zeros(3,length(nu));
    v = zeros(3,length(nu));

    % r and v for each true anomaly
    p = a*(1-e^2);
    for i = 1:length(nu)
        r(:,i) = p/(1+e*cos(nu(i))) * (cos(nu(i))*e_hat + sin(nu(i))*e_hat_perp);
        v(:,i) = sqrt(mu/p) * (-sin(nu(i))*e_hat + (e+cos(nu(i)))*e_hat_perp);
    end

end