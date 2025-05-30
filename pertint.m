%% ME512 Spaceflight Mechanics 
% Date: 10/12/2020
% Author: Nian Liu
% Description: variation in orbital elements of intermediate steps before 
% reaching final orbit due to perturbation around moon J2 & J3

function [ddt, d_f] = pertint(mu,e,a,om,inc,dt_hr)

% dt_hr = time in hour

R=1738E3; % radius of moon
J2 = 2.039E-4;% from class slides
J3 = 9.999E-6; 

%---circular for capture orbit---%
% mu = 4.905595800000000e+12;
% e = 0;
% a = 4536.09E3;
% om = 0;
% inc = deg2rad(63.43);
% dt = 0.2111 * 3600;
%---circular for capture orbit---%

dt = dt_hr * 3600; % convert time to seconds

n = sqrt(mu/a^3);

dedt = (3*n/2)*(R^3)/((a^3)*(1-e^2)^2) * J3 * sin(inc)*((5/4)*(sin(inc))^2-1) * cos(om);
didt = -(3*n/2) * (e*R^3)/((a^3)*(1-e^2)^3) * J3 * cos(inc)*((5/4)*(sin(inc))^2-1) * cos(om);
dOmdt = -(3*n/2) * cos(inc)* J2 * (R^2)/(a^2) - (3*n/2) * (e*R^3)/((a^3)*(1-e^2)^3) * J3 * cot(inc)*((15/4)*(sin(inc))^2-1) * sin(om);
domdt = (3*n/2) * J2 * (R^2) * (2 - 2.5* (sin(inc)^2)) /((a^2)*(1-e^2)) - ((3*n/2)*((J3 * R^3 * sin(om))/(2*e*a^3 * (1-e^2)^3 * sin(inc))))*(((5/4)*(sin(inc))^2-1)*(sin(inc))^2+e^2 * (1-(35/4)*(sin(inc))^2*(cos(inc))^2));

ddt = [dedt, didt, dOmdt, domdt]; % variation wrt time

de = dt * dedt;
di = dt * didt;
dOm = dt * dOmdt;
dom = dt * domdt;

d_f = [de, di, dOm, dom]; % variation in elements
