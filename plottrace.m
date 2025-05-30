%% ME512 Spaceflight Mechanics 
% Date: 08/12/2020
% Author: Nian Liu
% Description: Plot the trace of orbits in 3D (Earth reference frame cartesian) based on Keplerian orbital elements

function plottrace(kep_elements,mu,color_line,line_wid)

% kep_elements = 6 orbital elements [a e i Om om theta]
% mu = m_e; %commentout
% color_line = 'g-'; % commentout
% line_wid = 1.2; % commentout

for ii = 1:361
    
    kep_elements_loop = [kep_elements(1:5) deg2rad(ii)]; 

    kep_plot_all(ii,:) = kep2cart(kep_elements_loop,mu);

end

kep_plot = kep_plot_all(:,1:3);

plot3(kep_plot(:,1)/1000,kep_plot(:,2)/1000,kep_plot(:,3)/1000,color_line,'LineWidth',line_wid) % plot in km

grid on
axis equal