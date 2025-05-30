%% ME512 Spaceflight Mechanics 
% Date: 08/12/2020
% Author: Nian Liu
% Description: Plot the trace of orbits in 3D (translate from moon
% reference frame to Earth cartesian coordinates) from Keplerian orbital
% elements

function plottrace_m(kep_elements,mu,rotmat,r2,color_line,line_wid)

% kep_elements = 6 orbital elements [a e i Om om theta]
% r2; % position vector of moon
% mu = m_e; %commentout
% color_line = 'g-'; % commentout
% line_wid = 1.2; % commentout

for ii = 1:361
    
    kep_elements_loop = [kep_elements(1:5) deg2rad(ii)]; 

    kep_plot_all = kep2cart(kep_elements_loop,mu);
    
    kep_plot_all_m (ii,:) = rotmat\kep_plot_all(1:3).';% transform into earth ref frame

end

kep_plot = kep_plot_all_m(:,1:3);

plot3(r2(1)+kep_plot(:,1)/1000,r2(2)+kep_plot(:,2)/1000,r2(3)+kep_plot(:,3)/1000,color_line,'LineWidth',line_wid) % plot in km

grid on
axis equal