%% ME512 Spaceflight Mechanics 
% Date: 08/12/2020
% Author: Nian Liu
% Description: Find semimajor, right acension of ascending node, argument
% of pericenter and true anomaly based on statevectors (position, velocity)

function out = findascendingnode(r2, XVf, mu) % input unit m

% function [semimajor,lineofnode,ascendingnode,eccentricity,periapsisAnomaly,trueAnomaly] = findascendingnode(r2, XVf, mu)

% input in m

% semimajor

% mu = mu_e;
% r2 = r2*1000;
% XVf = 1000*XVf;

semimajor = mu*norm(r2)/(2*mu-norm(r2)*(norm(XVf))^2); %m

% line of nodes

h_angmom = crossproduct(r2, XVf);  % cross vector beween 2 vectors

k_polaraxis = [0 0 1];
j_axis = [0 1 0];
i_axis = [1 0 0];

lineofnode = crossproduct(h_angmom,k_polaraxis);

% argument of ascending node

if dotproduct(lineofnode,j_axis) > 0
    ascendingnode = acos(dotproduct(i_axis,lineofnode)/norm(lineofnode));
else
    ascendingnode = 2*pi() - acos(dotproduct(i_axis,lineofnode)/norm(lineofnode));    
end

% eccentricity

eccentricity_vec = (norm(XVf)^2*(r2) - (dotproduct(r2,XVf))*XVf)/mu - r2/norm(r2);
eccentricity = norm(eccentricity_vec);

% periapsis anomaly

if dotproduct(eccentricity_vec,k_polaraxis) >= 0
    periapsisAnomaly = acos(dotproduct(eccentricity_vec,lineofnode)/(eccentricity*norm(lineofnode)));
else
    periapsisAnomaly = 2*pi() - acos(dotproduct(eccentricity_vec,lineofnode)/(eccentricity*norm(lineofnode)));    
end

% true anomaly

semilatus = semimajor*(1-eccentricity^2);
radial_vel = dotproduct(XVf,r2)/norm(r2);

trueAnomaly_cos = semimajor*(1-eccentricity^2)/(eccentricity*norm(r2))-1/eccentricity;
trueAnomaly_sin = (radial_vel*semilatus)/(eccentricity*sqrt(semilatus*mu));

if trueAnomaly_sin >= 0
    trueAnomaly = acos(trueAnomaly_cos);
else
    trueAnomaly = 2*pi() - acos(trueAnomaly_cos);
end

out = [semimajor,ascendingnode,eccentricity,periapsisAnomaly,trueAnomaly];

