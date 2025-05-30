%% ME512 Spaceflight Mechanics 
% Date: 08/12/2020
% Author: Nian Liu
% Description: find inclination based on state vectors (position, velocity)

function i_actual = findinclination(r2, XVf)

h_angmom = crossproduct(r2, XVf); % use script for cross product
%h_angmom = [r2(2)*XVf(3)-r2(3)*XVf(2), r2(1)*XVf(3)-r2(3)*XVf(1), r2(1)*XVf(2)-r2(2)*XVf(1)];  % cross vector beween 2 vectors

k_polaraxis = [0 0 1];

i_cos = (h_angmom(1)*k_polaraxis(1)+h_angmom(2)*k_polaraxis(2)+h_angmom(3)*k_polaraxis(3))/norm(h_angmom);

% i_actual = rad2deg(acos(i_cos));

i_actual = acos(i_cos);

