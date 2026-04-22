function [r_E] = latlon2ECEF(r, phi, lambda)
%LATLON2ECEF Summary of this function goes here
%   Detailed explanation goes here

x = r*cosd(phi)*cosd(lambda);
y = r*cosd(phi)*sind(lambda);
z = r*sind(phi);

r_E = [x; y; z];
end
