%% ME512 Spaceflight Mechanics 
% Date: 08/12/2020
% Author: Nian Liu
% Description: calculates true anomaly travelled from perigee after given
% time by solving time equation numerically

function [trueAnomaly] = time2angle(time, mu, a, e)

% enter time in minutes
% a = semimarjor in meter

% %------------comment out
% time = 8000/60; % in min
% mu = mu_e;
% a = 14000E3; % example problem
% e = 0.4; % example problme
% %------------comment out

t = time * 60; % in sec
n = sqrt(mu/a^3);
M = n*(t);
E0 = M;
delta_E = (M - E0 + e*sin(E0))/(1 - e*cos(E0));

while norm(delta_E) > 0.01
   
    E0 = E0 + delta_E;
    delta_E = (M - E0 + e*sin(E0))/(1 - e*cos(E0));
    
end

trueAnomaly = 2 * atan(sqrt((1+e)/(1-e)) * tan(E0/2)); % in rad