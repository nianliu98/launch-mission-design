%% ME512 Spaceflight Mechanics 
% Date: 11/12/2020
% Author: Nian Liu
% Description: Main script that calculates orbit transfer relevant information
% Other function files required: cart2kep, crossproduct, dotproduct,
% findascendingnode, findinclination, kep2cart, lambertMR, moon_eph,
% plottrace, plottrace_m, time2angle, pertint

clc
clear
close all

%------Changeable parameters in this code------%

%%% line 59: if P.O. is leo, leo altitude
%%% line 84: date/time of arrival at moon
%%% line 108-112: SPS moon orbit parameters
%%% line 199: flight angle for transfer
%%% line 204: Delayed time in gto for delayed prop mass calc
%%% line 208: Parking orbit can be gto/leo by changing semimarjor variable
%   between gto_a and r_1
%%% line 390: delta V trim
%%% line 453-459: Uncomment for prop mass for full year, line 83 j = 364, set line 84 MJD = [2030,1,1,12,0,0]
%%% line 470-475: Uncomment for porkchop, line 82-83 loop j values need > 1

%%%------Changeable parameters in this code------%%%


%% Earth and Moon Parameters

% _m: moon, _e: earth

% mass
m_m = 7.35E22;
m_e = 5.972E24;

% G
G = 6.67428E-11;
g_0 = 9.81; 

% gravitational constant mu
mu_m = G * m_m;
mu_e = G * m_e;

% earth radius (m)
r_e = 6378E3;

% moon orbit info
a_m = 384400E3; % m, semi-major axis
r_m = 1738E3; % m, radius of moon
peri_m = 363300E3; %perigee
apo_m = 405500E3; %apogee
e_m = 0.0549; % moon eccentricity
i_m = 28.54; %moon orbit plane inclination relative to earth equatorial



%% Input Parameters
parkingOrbitAlt = 1000E3;

    % calculate delta v1 required injecting from LEO to Hohmann transfer ellipse perigee
    r_1 = r_e + parkingOrbitAlt; %leo radius
    v_leo = sqrt(mu_e/r_1); %velocity at circular leo with 0 inclination
    v_htp = sqrt((2*mu_e/r_1)-(2*mu_e/(r_1+peri_m))); %velocity at Hohmann transfer ellipse perigee

    delta_v1 = v_htp - v_leo;
    
    % calculate delta v2 required injecting from Hohmann transfer ellipse apogee to moon orbit with change in inclination
    r_2 = peri_m; %assuming injection to moon orbit perigee
    v_m = sqrt((2*mu_e/r_2)-(mu_e/a_m)); %calculating moon velocity at moon orbit perigee
    v_hta = sqrt((2*mu_e/r_2)-(2*mu_e/(r_1+peri_m))); %velocity at Hohmann transfer ellipse apogee

    delta_v2 = sqrt(v_m^2 + v_hta^2 - 2 * v_m * v_hta * cos(deg2rad(i_m)));

    
    %% moon capture
    
    %%% doing the capture asysmptote thing

% final modified julian dates (moon perigee?)

for jj = 0
for j = 0
mjd2000f = mjuliandate([2030,2,26,12,0,0])+j-51544.5; % +98 is 99th day

% final moon position & velocity (where we intersect, perigee)
[r2,XVf] = moon_eph(mjd2000f); % km

% get moon orbital elements
kep_ele_m = cart2kep(r2*1000, XVf*1000, mu_e); % [a e i Om om theta]

% [semimajor,lineofnode,angleoffset,eccentricity,periapsisAnomaly,trueAnomaly] = findascendingnode(r2*1000, XVf*1000, mu_e); % get all orbit elements
% i_mcalc = findinclination(r2, XVf);
semimajor = kep_ele_m(1);
i_mcalc = kep_ele_m(3);
angleoffset = kep_ele_m(4);
eccentricity = kep_ele_m(2);
periapsisAnomaly = kep_ele_m(5);
trueAnomaly = kep_ele_m(6);

    refangle = deg2rad(23.4-1.54); % polar axis of moon ref rotation relative to earth polar axis
    rotmat = [1 0 0; 0 cos(refangle) sin(refangle); 0 -sin(refangle) cos(refangle)]; % rotation matrix

% period of moon
period_m = 2*pi()*sqrt(semimajor^3/mu_e);
    
%%% SPS orbital parameters
    alt_sps = 2798.09E3; %196.3E3; %400E3; not too much difference in inclination change delta v for diff altitude
    e_orb = 0.3;
    r_sps = r_m + alt_sps; %m, radius of sps orbit around moon from moon centre
    a_sps = r_sps; % assuming circular orbit, semimajor axis of sps orbit
    inc_sps = 63.43;
    
    %%% SPS operational plane circular orbit (first captured) period
    T_sps_cap = 2*pi()*sqrt(r_sps^3/mu_e)/3600; % hr
     
    %%% time for sps from capture at 90 deg to 180 deg for peri lowering
    t_sps_90to180 = T_sps_cap/4;
    disp('Time for SPS to travel from 90 capture to 180 peri-lowering (hr):')
    disp(t_sps_90to180)
    
    %%% time for lander from inc change at TA 90 deg to TA 360 deg circ orbit to
    %%% decent
    disp('Time for lander to travel from 90 capture to 360 decent (hr):')
    disp(t_sps_90to180*3)

    % sps orbit when first entering moons sphere of influence around moon
    % in plane with orbit transfer)
    
    % capture orbit argument of pericentre om
    om_capture = 0;
%     om_capture = periapsisAnomaly+trueAnomaly+pi(); %rad
    theta_capture = 10; %deg
    
    if inc_sps > 75
        i_capture = 70; % deg. inclination of capture orbit around moon, moon ref
    else
        i_capture = inc_sps;
    end
    
    i_sps = inc_sps - i_capture; % sps operational orbit relative to cpature orbit plane

    kep_sps_equ_local = ([r_sps 0 deg2rad(i_capture) angleoffset om_capture deg2rad(theta_capture)]); %m, [a e i Om om theta]
%     kep_sps_equ_local = ([r_sps 0 deg2rad(i_capture) deg2rad(180) deg2rad(om_capture) deg2rad(theta_capture)]); %m, [a e i Om om theta]
    kep_sps_equ = ([r_sps 0 deg2rad(i_capture) deg2rad(180) om_capture deg2rad(theta_capture)]); %m, [a e i Om om theta]

    state_cap_m=kep2cart(kep_sps_equ_local,mu_m); % m, moon ref frame, sps state around moon capture orbit
    rs_cap_m = state_cap_m(1:3);
    rs_cap_e = transpose(rotmat\rs_cap_m.'); %m, transform to earth ref frame
     % position vector of sps capture relativ to moon in earth ref
    vs_cap_m = state_cap_m(4:6); % velocity vector of sps capture relativ to moon in moon ref
    vs_cap_e = transpose(rotmat\vs_cap_m.'); % velocity vector of sps capture relativ to moon in earth ref
    
    rf = r2*1000 + rs_cap_e; % m, position vector of sps capture relative to Earth in earth ref
   
    
    [x,y,z] = sphere;
    xs_cap = x * 800;
    ys_cap = y * 800;
    zs_cap = z * 800;
    surf(xs_cap+rf(1)/1000,ys_cap+rf(2)/1000,zs_cap+rf(3)/1000,'FaceAlpha',0.5) %km
    
    hold on
    
    v_h = sqrt(2*mu_m/r_sps+delta_v2^2); %hyperbolic velocity
    v_sps = sqrt(2*mu_m/r_sps-mu_m/a_sps); % sps velocity around moon
    
    % moon orbit plane change to enter operational orbit
    % velocity stays the same before and after plane change
    
    %  performs plane change during capture to enter operational orbit
    delta_c = sqrt(v_h^2 + v_sps^2 - 2 * v_h * v_sps * cos(deg2rad(i_sps)));
    
    
    % total delta v budget
    delta_vtot = delta_v1 + delta_c;  
    
%     disp('delta V capture')
%     disp(delta_c)
%     
%     disp('Total deltaV budgetfor Hohmann is')
%     disp(delta_vtot)
    

%% Lambert


%%% earth stuff based on moon orbit
%%% define gto orbit

gto_apo = r_e + 35786E3; %m, apogee radius
gto_peri = r_e + 185E3; %m, perigee radius
gto_a = (gto_apo + gto_peri)/2; % gto semiajor axis
gto_e = (gto_apo - gto_peri)/(gto_apo + gto_peri); % gto eccentricity

periapsis_lat = periapsisAnomaly + trueAnomaly; % argument of periapsis lattitude

% find a good theta angle
flightangle = 178; %deg angle beteween moon vector and sat orbit

periapsis_lat_om = periapsis_lat + deg2rad(360 - flightangle); % offet the gto perigee by a 
% bit so when flight angle is not 180, still transferring from gto apogee

delaytime = 0; % min
trueAnomaly_delay = time2angle(delaytime, mu_e, gto_a, gto_e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kep_leo = ([gto_a gto_e i_mcalc angleoffset periapsis_lat_om 0]); %m, [a e i Om om theta] for gto
kep_leo = ([gto_a gto_e i_mcalc angleoffset periapsis_lat_om trueAnomaly_delay]); %m, [a e i Om om theta] for gto delayed
% kep_leo = ([r_1 0 i_mcalc angleoffset periapsis_lat_om 0]); %m, [a e i Om om theta] for leo

% initial conditions
stateinitial=kep2cart(kep_leo,mu_e); % m, output contains State vector in cartesian coordinates (position,velocity)

    % Time takes for Hohmann transfer
    a_h = (gto_peri+norm(r2)*1000)/2; % semi-major axis of Hohmann transfer ellipse
    t_h = pi()*sqrt(a_h^3/mu_e); % half of Homann ellipse period (s)
    t_hhr = t_h/3600; % time in hour
    disp('Time takes for Hohmann transfer (hr)')
    disp(t_hhr)

    % time of flight
    tofh = 80 + 4 * 5;
%     tofh = 80 + 4 * 5 - delaytime/60 ; %hour, for delay on ground with nominal arrival
    tof = tofh * 3600; %s

% period of leo
period_leo = 2*pi()*sqrt(gto_a^3/mu_e); %s

% display earth (km)

xe = x * r_e./(10^3);
ye = y * r_e./(10^3);
ze = z * r_e./(10^3);

surf(xe,ye,ze,'FaceAlpha',0.5)
text(-10000,0,40000,'Earth')
axis equal


% display final position of moon 

xm = x * r_m./(10^3); % km
ym = y * r_m./(10^3);
zm = z * r_m./(10^3);

% xm = 2* x * r_m./(10^3); % km % moon scaled up by 2
% ym = 2* y * r_m./(10^3);
% zm = 2* z * r_m./(10^3);

% % plot original scale
% surf(xm+r2(1),ym+r2(2),zm+r2(3),'FaceAlpha',0.5) % Earth moon distance scaled by factor of 10
% plot3(xout,yout,zout,'b-','LineWidth',2) % plot leo orbit
% plot3(xout_m,yout_m,zout_m,'r-','LineWidth',1) % plot moon orbit scale by factor of 10

% scale down by factor of 10
% surf(xm+r2(1)/10,ym+r2(2)/10,zm+r2(3)/10,'FaceAlpha',0.5) % Earth moon distance scaled by factor of 10
% text(r2(1)/10,r2(2)/10,r2(3)/10+5000,'Moon')

surf(xm+r2(1),ym+r2(2),zm+r2(3),'FaceAlpha',0.5) % Earth moon distance not scaled by factor of 10
text(r2(1),r2(2),r2(3)+30000,'Moon')


% display satellite initial position
% make sat sphere
xs = x * 800; % km
ys = y * 800;
zs = z * 800;

surf(xs+stateinitial(1)/1000,ys+stateinitial(2)/1000,zs+stateinitial(3)/1000,'FaceAlpha',0.5) %km

% plot orbits

plottrace(kep_leo,mu_e,'b-',1.2) % plot leo orbit
% plot3(xout_m/10,yout_m/10,zout_m/10,'r-','LineWidth',1) % plot moon orbit scale by factor of 10
plottrace(kep_ele_m,mu_e,'r-',1.2) % plot moon orbit
plottrace_m(kep_sps_equ_local,mu_m,rotmat,r2,'m-',1.2) % plot sps capture orbit

xlabel('x(km)')
ylabel('y(km)')
zlabel('z(km)')

% use lambert
r1 = stateinitial(1:3); % m, position vector of module in parking orbit
v1 = stateinitial(4:6); % m/s, velocity vector of module in parking orbit

v1_mag = norm(v1);

[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r1,r2*1000,tof,mu_e,0,0,0); % meet moon centre

% kep_ele_tr = cart2kep(r2*1000,VF, mu_e); % meet moon centre
kep_ele_tr = cart2kep(r2*1000,VF, mu_e); % meet moon SoI

%%% plot transfer orbit
kep_transfer = ([A E kep_ele_tr(3) kep_ele_tr(4) kep_ele_tr(5) deg2rad(0)]); %m, [a e i Om om theta]
% kep_transfer = ([A E deg2rad(i_mcalc) deg2rad(angle_lon+180) deg2rad(282.2) deg2rad(0)]); %m, [a e i Om om theta]

% period of transfer
period_transfer = 2*pi()*sqrt(A^3/mu_e); %s

%------------------------------------------------
plottrace(kep_transfer,mu_e,'g-',1) % plot transfer orbit

v2 = XVf*1000; % m/s, velocity vector of moon

DELTAV1 = VI - v1;
DELTAV2 = v2 - VF; % use moon vel
% DELTAV2 = vs_cap_e - VF; %v2-VF for intersecting moon instead of SoI
DELTAV1_mag = norm(DELTAV1);
DELTAV2_mag = norm(DELTAV2);

%%% declination of asymptote
% change v infinity (delta v2) into moon ref frame
DELTAV2_mref = rotmat * transpose(DELTAV2);
DELTAV2_mref_mag = norm(DELTAV2_mref);
declination = asin(norm(DELTAV2_mref(3))/DELTAV2_mref_mag); % rad, v infnitiy in moon ref frame

i_hyperbola = rad2deg(asin(sin(declination)/(sin(deg2rad(rad2deg(om_capture) + theta_capture))))); % deg, inclination of hyperbola plane
vsps_cap = state_cap_m(4:6); % moon ref

DELTAV_HYPER = sqrt(2*mu_m/r_sps + DELTAV2_mref_mag^2);

% use DELTAV2_mref as v infinity, state of sps orbit around moon capture as
% v capture, the delta v required while performing an inclination hyperbola
% change is:
DELTAV_CAP = sqrt(DELTAV_HYPER^2 + norm(vsps_cap)^2 - 2 * DELTAV_HYPER * norm(vsps_cap) * cos(deg2rad(i_capture-i_hyperbola)));

%%% total delta v
% delta v for sending sps and base together orbit transfer into capture
% orbit inclination
DELTAV_tot_spsbase = norm(DELTAV1)+DELTAV_CAP; % transfer delta V for sps + base
disp('Total transfer delta V:')
disp(DELTAV_tot_spsbase)

% lunar base module plane change, from the capture orbit plane,
% i_sps is the inclination change required for sps, but module has a lower
% inclintation, so subtracted from sps inc change required

vbase_cap_mag = norm(vsps_cap);
% dv for sps + base into landing orb if inc_sps > 75, 75 Shrodinger 75S
% dv for just base into landing orb if inc_sps < 75
DELTAV_inc_base = vbase_cap_mag*sin(deg2rad(75 - i_capture)/2); 

% SPS plane change into polar orbit from orbit of lunar base, need further inclination
% change of 15 deg, difference between polar 90 and lunar base orbit 75

if inc_sps > 75 % if sps operational orbit > 75S, need another inc change to enter operational
    
    DELTAV_inc_sps = 2*vbase_cap_mag*sin(deg2rad(inc_sps-75)/2); % dv for sps to enter operational orbit after lander landed

else
    
    DELTAV_inc_sps = 0; % sps already captured at operational orbit, no need for another inc change
    
end
%% Landing maneuvre
%%% take the procedure similar to hohmann transfer

I_sp = 340; %PPS-1350

a_land = (r_sps+r_m)/2;
e_land = (r_sps-r_m)/(r_sps+r_m);
deltav1_land = sqrt(2*mu_m/r_m-mu_m/a_land); % for bringing lander to rest on surface of moon
deltav2_land = sqrt(mu_m/r_sps) - sqrt(2*mu_m/r_sps-mu_m/a_land); %required for entering the decsending orbit
deltav_tot_land = deltav1_land + deltav2_land;
m_f_base = 2000; % dry mass of lander
m_i_base = m_f_base*(exp((deltav_tot_land)/(I_sp*g_0)));
m_p_land = m_i_base - m_f_base; % propellant required to land module

% time taken for landing maneuvre (approx hohmann)
T_landing = pi()*sqrt((a_land^3)/mu_m); % in s
T_landing_h = T_landing/3600; % in hour
disp('Time taken for landing maneuvre (hr):')
disp(T_landing_h)

%% propellant


% delta v for sps entering elliptical op orbit from apogee

ra_orb = r_sps;
rp_orb = ra_orb * (1 - e_orb);
a_orb = (rp_orb + ra_orb)/2;

v_orb_ini = sqrt(mu_m/ra_orb);
v_orb_a = sqrt(2*mu_m/ra_orb - mu_m/a_orb);

DELTAV_sps_ellip = norm(v_orb_a-v_orb_ini);

% orbit trim

delta_v_trim = 130.2925; % for 5 years
m_f_trim = 7000;
m_i_trim = m_f_trim*(exp((delta_v_trim)/(I_sp*g_0)));
m_p_trim = m_i_trim - m_f_trim;

%%% prop mass for sps entering op orbit
m_f_sps = 7000 + m_p_trim;
m_i_sps = m_f_sps*(exp((DELTAV_inc_sps+DELTAV_sps_ellip)/(I_sp*g_0))); % propellant for sps to enter operational orbit
m_p_sps = m_i_sps - m_f_sps;


%%% prop mass for just lander entering landing 75S orbit if in_sps < 75
if inc_sps < 75
    m_f_landorbit = 2000 + m_p_land;
    m_i_landorbit = m_f_landorbit*(exp((DELTAV_inc_base)/(I_sp*g_0))); % propellant for sps to enter operational orbit
    m_p_landorbit = m_i_landorbit - m_f_landorbit;
else
    m_p_landorbit = 0;
end

%%% prop mass for sps + lander
m_f_spsbase = 9000 + m_p_land + m_p_sps + m_p_landorbit + m_p_trim; % mass carried including propellant for landing & enter sps orbit

if inc_sps > 75
    m_i_spsbase = m_f_spsbase*(exp((DELTAV_tot_spsbase+DELTAV_inc_base)/(I_sp*g_0))); % prop for sps+bas transfer & capture change to 75S
else
    m_i_spsbase = m_f_spsbase*(exp((DELTAV_tot_spsbase)/(I_sp*g_0))); % prop for base change to 75S
end

m_p_spsbase = m_i_spsbase - m_f_spsbase; % propellent required for sps and base getting into landing orbit


%% total propellant required excluding orbit trim

% for changing plane in leo

if i_mcalc < 28.5
    DELTAV_leo = 2*v_leo*sin(deg2rad((28.5-i_mcalc)/2));
    m_f_everything = m_p_spsbase + m_p_sps + m_p_land + m_p_trim + 9000;
    m_i_everything = m_f_everything *(exp((DELTAV_leo)/(I_sp*g_0)));
    m_p_everything = m_i_everything - m_f_everything;
else
    m_p_everything = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m_tot_prop(jj+1,j+1) = m_p_spsbase + m_p_sps + m_p_land + m_p_landorbit + m_p_trim + m_p_everything; %kg, for leo
m_tot_prop(jj+1,j+1) = m_p_spsbase + m_p_sps + m_p_land + m_p_landorbit + m_p_trim ; %kg, for gto

disp('Total propellant required:')
disp(m_tot_prop)

disp('Total Launch mass:')
disp(m_tot_prop+9000)

% i_mcalc_everyday(j+1) = i_mcalc;

hold off

end
tofh_plot(jj+1) = tofh;
end

%----- plot total prop mass in year 2030 -----%
% figure
% plot(1:j+1,m_tot_prop)
% grid on
% xlabel('Day of Year 2030')
% ylabel('Total Required Propellant Mass (kg)')
%----- plot total prop mass in year 2030 -----%

min_prop = min(min(m_tot_prop))
[row_prop, column_prop] = find(m_tot_prop == min_prop); % find day of year with min propellant mass

disp('day of lowest prop mass')
disp(column_prop)

disp('True anomaly delayed:')
disp(rad2deg(trueAnomaly_delay))

%----- pork chop ----%
% figure
% contourf(10:j+10,tofh_plot,m_tot_prop,'ShowText','on')
% xlabel('Arrival Day in 2030')
% ylabel('Time of Flight (h)')
%----- pork chop ----%

%% Delay in orbit

% if delayed in parking orbit

% change true anomaly of injection in previous scripts

%% Parking orbit secular variation

J2 = 0.0010826266;
n_pert = sqrt(mu_e/(gto_a^3));
j2 = (3/2) * J2 * (r_e/(gto_a*(1-gto_e^2)))^2;

nj2 = n_pert*j2;
% nj2 = 0.96404 * (r_e/r_1)^(7/2); % deg/day, for cirbutlar??

% regression of line of nodes
Om_sec_rate = -nj2 * cos(i_mcalc); % deg, change over one period
om_sec_rate = 2*nj2*(1-(5/4)*(sin(i_mcalc))^2);

% over one period, delta v required to correct small omega
% dont need to correct argument of periapsis because circular orbit


%% Parking orbit drag

% payload diameter of Falcon Heavy
d_payload = 5.2; %m
% cross-sectional area
area_cross = pi()*(d_payload/2)^2;
% total wet mass
m_tot_wet = m_tot_prop + 9000;
% assumed coefficient
C_d = 2.2;
% ballistic coefficient
BC = (C_d*area_cross)./m_tot_wet;
% BC = (2*1529.6)/419700; %iss
%---------
% scaled Height of Earth's atomosphere
H_scaled = (1.38 * 10^(-23) * 290)/(g_0 * 4.76*10^(-26)); %m
% air density at desired altitude
rho_alt = 1.225*exp(-(185E3/H_scaled)); % at gto perigee
% mean motion
n_meanmotion = sqrt(mu_e/gto_a^3);
% rate of change of semi-major axis, hence radius
sigma_e_theta = sqrt(mu_e/gto_a)* ( sqrt(1-gto_e^2) / v1_mag);
a_dot = -(BC*rho_alt*v1_mag^2)/(n_meanmotion * sigma_e_theta);

% over one period
a_dot_period = a_dot * period_leo;
deltav_drag = (n_meanmotion/4)*sqrt((1+gto_e)/(1-gto_e)).*a_dot;

% a_t = (sqrt(r_1) - sqrt(mu_e)*rho_alt*BC*86400).^2 - r_1;

%% Perturbations in intermediate orbits

% d_f are final variations in [e, i, Om, om]

%%% SPS 0.211hr in capture plane
[ddt1, d_f1] = pertint(mu_m,0,r_sps,0,deg2rad(i_capture),t_sps_90to180)
%%% lander 0.633 hr in landing circular orbit
[ddt2, d_f2] = pertint(mu_m,0,r_sps,0,deg2rad(75),t_sps_90to180*3)
%%% lander 2.189 hr in landing elliptical orbit
[ddt3, d_f3] = pertint(mu_m,e_land,a_land,deg2rad(180),deg2rad(75),T_landing_h)


%%
% ax = gca;               % get the current axis
% ax.Clipping = 'off';    % turn clipping off