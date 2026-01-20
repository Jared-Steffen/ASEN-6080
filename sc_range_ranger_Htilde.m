function [Htilde] = sc_range_ranger_Htilde(sc_state,gs_state)
%{
Inputs:
    >sc_state: 6x1 S/C state vector
    >gs_state: 6x1 GS state vector
Outputs:
    >Htilde: jacobian of partials of rho and rho_dot wrt S/C state
%}

% Break apart state vectors
x = sc_state(1);
y = sc_state(2);
z = sc_state(3);
xd = sc_state(4);
yd = sc_state(5);
zd = sc_state(6);
xs = gs_state(1);
ys = gs_state(2);
zs = gs_state(3);
xsd = gs_state(4);
ysd = gs_state(5);
zsd = gs_state(6);

% Assign differnce vars
dx = x-xs;
dy = y-ys;
dz = z-zs;
dxd = xd - xsd;
dyd = yd - ysd;
dzd = zd - zsd;
rho = norm(sc_state(1:3)-gs_state(1:3));
N = dx*dxd + dy*dyd + dz*dzd;

Htilde = [dx/rho, dy/rho, dz/rho,...
          0, 0, 0;
          dxd/rho-(dx*N)/rho^3, dyd/rho-(dy*N)/rho^3, dzd/rho-(dz*N)/rho^3,...
          dx/rho, dy/rho, dz/rho];

end

