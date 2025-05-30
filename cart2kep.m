%% ME512 Spaceflight Mechanics 
% Date: 08/12/2020
% Author: Nian Liu
% Description: Convert statevectors (position and velocity) into Kepler
% orbital elements

function kep_ele = cart2kep(r2, XVf, mu)

% ascend_ele = [semimajor,ascendingnode,eccentricity,periapsisAnomaly,trueAnomaly]

ascend_ele = findascendingnode(r2, XVf, mu);

% find elements

kep_a = ascend_ele(1);
kep_e = ascend_ele(3);
kep_i = findinclination(r2, XVf);
kep_Om = ascend_ele(2);
kep_om = ascend_ele(4);
kep_theta = ascend_ele(5);

kep_ele = [kep_a kep_e kep_i kep_Om kep_om kep_theta]; %[a e i Om om theta]