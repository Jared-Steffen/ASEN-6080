function [phi] = STM_J2_J3(mu,J2,J3,R,r_vec)
%{
NOTE: This considers the effects of J2 and J3 orbital pertubations
Inputs:
    >mu: gravitational parameter of central body (km^3/s^2)
    >J2: J2 coefficient of central body
    >J3: J3 coefficient of central body
    >R: radius of central body
    >r_vec: 3D cartesian S/C position relative to central body
Outputs:
    >stm: STM made of partials of acceleration wrt S/C position
          and mu, J2, J3 (9x9)
%}

% Break apart r vector and get magnitude
x = r_vec(1);
y = r_vec(2);
z = r_vec(3);
r = norm(r_vec);

% Initialize
phi = zeros(9);

% Fill in dv/dv portion
phi(1:3,4:6) = eye(3);

% Fill in grav coefficients section
phi(4:6,7) = -1/r^3 .* [x, y, z]' + ...
    (3*J2*R^2)/(2*r^5) .* [x*(5*z^2/r^2-1), y*(5*z^2/r^2-1), z*(5*z^2/r^2-3)]'...
    + (J3*R^3)/(2*r^5) .* [5*x/r*(7*z^3/r^3-3*z/r), 5*y/r*(7*z^3/r^3-3*z/r), 35*z^4/r^4-30*z^2/r^2+3]'; % mu

phi(4:6,8) = (3*mu*R^2)/(2*r^5) .* [x*(5*z^2/r^2-1), y*(5*z^2/r^2-1), z*(5*z^2/r^2-3)]'; % J2

phi(4:6,9) = (mu*R^3)/(2*r^5) .* [5*x/r*(7*z^3/r^3-3*z/r), 5*y/r*(7*z^3/r^3-3*z/r), 35*z^4/r^4-30*z^2/r^2+3]'; % J3

% Acceleration partials

da_mudr = mu/r^5.*[-(r^2-3*x^2), 3*x*y, 3*x*z;
                   3*x*y, -(r^2-3*y^2), 3*y*z;
                   3*x*z, 3*y*z, -(r^2-3*z^2)];

da_J2dr = [(3*J2*mu*R^2)/(2*r^5)*(5*z^2/r^2*(1-7*x^2/r^2)-(1-5*x^2/r^2)),...
           (15*J2*mu*R^2*x*y)/(2*r^7)*(1-7*z^2/r^2),...
           (15*J2*mu*R^2*x*z)/(2*r^7)*(3-7*z^2/r^2);
           (15*J2*mu*R^2*x*y)/(2*r^7)*(1-7*z^2/r^2),...
           (3*J2*mu*R^2)/(2*r^5)*(5*z^2/r^2*(1-7*y^2/r^2)-(1-5*y^2/r^2)),...
           (15*J2*mu*R^2*y*z)/(2*r^7)*(3-7*z^2/r^2);
           (15*J2*mu*R^2*x*z)/(2*r^7)*(3-7*z^2/r^2),...
           (15*J2*mu*R^2*y*z)/(2*r^7)*(3-7*z^2/r^2),...
           (3*J2*mu*R^2)/(2*r^5)*(5*z^2/r^2*(3-7*z^2/r^2)-3*(1-5*z^2/r^2))];

da_J3dr = [(5*J3*mu*R^3*z)/(2*r^7)*(7*z^2/r^2*(1-9*x^2/r^2)-3*(1-7*x^2/r^2)),...
           (105*J3*mu*R^3*x*y*z)/(2*r^9)*(1-3*z^2/r^2),...
           (5*J3*mu*R^3*x)/(2*r^7)*(7*z^2/r^2*(3-9*z^2/r^2)-3*(1-7*z^2/r^2));
           (105*J3*mu*R^3*x*y*z)/(2*r^9)*(1-3*z^2/r^2),...
           (5*J3*mu*R^3*z)/(2*r^7)*(7*z^2/r^2*(1-9*y^2/r^2)-3*(1-7*y^2/r^2)),...
           (5*J3*mu*R^3*y)/(2*r^7)*(7*z^2/r^2*(3-9*z^2/r^2)-3*(1-7*z^2/r^2));
           -(15*J3*mu*R^3*x)/(2*r^7)*(21*z^4/r^4-14*z^2/r^2+1),...
           -(15*J3*mu*R^3*y)/(2*r^7)*(21*z^4/r^4-14*z^2/r^2+1),...
           (5*J3*mu*R^3*z)/(2*r^7)*(7*z^2/r^2*(4-9*z^2/r^2)-6*(2-7*z^2/r^2)-3)];

phi(4:6,1:3) = da_mudr + da_J2dr + da_J3dr;

end


