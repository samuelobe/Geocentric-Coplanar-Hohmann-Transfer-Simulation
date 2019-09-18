% Name: Samuel Obe
% Project: Geocentric Coplanar Hohmann Transfer Simulation

r_earth = 6378.137; % Equatorial radius of the Earth (km)
mu_earth = 398600; % Standard gravitational parameter of the Earth (km^3/s^2)

fprintf('Geocentric Coplanar Hohmann Transfer Simulation\n');
fprintf('\nRequired inputs:\n');

while(1)
    circ_1_alt = input('\nPlease input the altitude for the initial circular orbit (km): \n');
    if circ_1_alt < 160
        fprintf('ERROR: UNSTABLE ORBIT PLEASE TRY AGAIN!')
    else
        circ_1_r = circ_1_alt + r_earth;
        break
    end 
     
end

while(1)
    circ_2_alt = input('\nPlease input the altitude for the final circular orbit (km): \n');
    if (circ_2_alt+r_earth) <= circ_1_r
        fprintf('ERROR: FINAL CIRCULAR ORBIT IS LESS THAN INITIAL ORBIT PLEASE TRY AGAIN!')
    else
        circ_2_r = circ_2_alt + r_earth;
        break
    end 
     
end
%% Circular Orbit Velocities

circ_1_v = sqrt(mu_earth/circ_1_r); % Circular velocity of the initial orbit (km/s)
circ_2_v = sqrt(mu_earth/circ_2_r); % Circular velocity of the final orbit (km/s)

%% Transfer Orbit Calculations

T_eccent = (max(circ_1_r, circ_2_r) - min(circ_1_r, circ_2_r))/(circ_1_r + circ_2_r); % Transfer orbit eccentricity 
T_semi = (circ_1_r+circ_2_r)/2; % Transfer orbit semi major axis (km)
T_p_r = T_semi * (1-T_eccent); % Transfer orbit perigee distance (km)
T_a_r = T_semi * (1+T_eccent); % Transfer orbit apogee distance (km)
T_p_v = sqrt((2*mu_earth*circ_2_r)/(circ_1_r*(circ_1_r+circ_2_r))); % Transfer orbit perigee velocity (km/s)
T_a_v = sqrt((2*mu_earth*circ_1_r)/(circ_2_r*(circ_1_r+circ_2_r))); % Transfer orbit apogee velocity (km/s)
tof = (((2*pi/sqrt(mu_earth))*((T_p_r+T_a_r)/2)^(3/2))/2)/3600; % Time of flight (hours)

%% Delta V Calculations
delta_v_1 = T_p_v - circ_1_v;
delta_v_2 = circ_2_v - T_a_v;
delta_v_total = delta_v_1 + delta_v_2;

%% Outputs
fprintf('First Circular Orbit Altitude = %.4f km \n\n', circ_1_alt);
fprintf('First Circular Orbit Radius = %.4f km \n\n', circ_1_r);
fprintf('First Circular Orbit Velocity = %.4f km/s \n\n', circ_1_v);
fprintf('Second Circular Orbit Altitude = %.4f km \n\n', circ_2_alt);
fprintf('Second Ciruclar Orbit Radiius = %.4f km \n\n', circ_2_r);
fprintf('Second Ciruclar Orbit Velocity = %.4f km/s \n\n', circ_2_v);
fprintf('Transfer Orbit Eccentricity = %.4f \n\n', T_eccent);
fprintf('Transfer Orbit Semimajor Axis = %.4f km \n\n', T_semi);
fprintf('Transfer Orbit Perigee Distance = %.4f km \n\n', T_p_r);
fprintf('Transfer Orbit Apogee Distance = %.4f km \n\n', T_a_r);
fprintf('Transfer Orbit Velocity at Perigee = %.4f km/s \n\n', T_p_v);
fprintf('Transfer Orbit Velocity at Apogee = %.4f km/s \n\n', T_a_v);
fprintf('Time of Flight = %.4f hours \n\n', tof);
fprintf('Delta V 1 = %.4f km/s \n\n', delta_v_1);
fprintf('Delta V 2 = %.4f km/s \n\n', delta_v_2);
fprintf('Delta V total = %.4f km/s \n\n', delta_v_total);

%% Settings for 3D plot
figure(1)
hold on
grid on
axis equal
rotate3d on

%% Plot the Earth
h1 = gca; 
earth_sphere(h1,'km') 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

%% Plot Circular Orbits
number_of_points = 1000;
theta = linspace(0,2*pi,number_of_points);
rho1 = ones(1,number_of_points)*circ_1_r;
rho2 = ones(1,number_of_points)*circ_2_r;
[X1,Y1] = pol2cart(theta,rho1);
[X2,Y2] = pol2cart(theta,rho2);
Z1 = 0*ones(1,length(X1));
Z2 = 0*ones(1,length(X1));

plot3(X1,Y1,Z1, '-c', 'LineWidth', 1.5);
plot3(X2,Y2,Z2, '-m', 'LineWidth', 1.5);

%% Plot Transfer Orbit
oe_transfer(1) = T_semi;
oe_transfer(2) = T_eccent;
oe_transfer(3) = 0;
oe_transfer(4) = 0;
oe_transfer(5) = 0;
oe_transfer(6) = 0;

[r_transfer, v_transfer] = orb2eci(mu_earth, oe_transfer);

period = 2.0 * pi * oe_transfer(1) * sqrt(oe_transfer(1) / mu_earth);
simulation_time = -(0.5*period/360); 

for i = 1:1:361
    simulation_time = simulation_time + (0.5*period/360);
    [r, ~] = twobody2 (mu_earth, simulation_time, r_transfer, v_transfer);
    r_x(i) = r(1);
    r_y(i) = r(2);
    r_z(i) = r(3);

end

plot3(r_x, r_y, r_z, '-r', 'LineWidth', 1.5);
plot3(r_x(1), r_y(1), r_z(1),'ok','MarkerSize',7,'MarkerFaceColor','k')
plot3(r_x(end), r_y(end), r_z(end),'ok','MarkerSize',7,'MarkerFaceColor','k')


