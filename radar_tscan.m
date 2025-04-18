clear; close all; clc;
% Constants and Earth parameters
videoFile = VideoWriter('Radar_simulation.mp4', 'MPEG-4');
videoFile.FrameRate = 60; % Set frame rate (adjust as needed)
open(videoFile); % Open the video file for writing
set(gcf, 'Renderer', 'opengl'); % High-quality rendering

Re = 6378;                % Earth's radius in km
mu_earth = 398600;              % Earth's gravitational parameter [km^3/s^2]

% Altitude range [km]
altitude_min = 600;
altitude_max = 645;

% Inclination range (degrees)
inclination_min = -75;
inclination_max = 75;

% RAAN range (degrees)
azimuth_min = 0;
azimuth_max = 360;

% Convert angles to radians
inclination_min_rad = deg2rad(inclination_min);
inclination_max_rad = deg2rad(inclination_max);
azimuth_min_rad = deg2rad(azimuth_min);
azimuth_max_rad = deg2rad(azimuth_max);

% Define radii
radius_min = Re + altitude_min;
radius_max = Re + altitude_max;

% Create spherical coordinate grid
[inclination, azimuth] = meshgrid(linspace(inclination_min_rad, inclination_max_rad, 100), ...
                                   linspace(azimuth_min_rad, azimuth_max_rad, 100));



%% SATELLITE ORBIT
a = Re + 595;  % Semi-major axis (km)
e = 0;              % Eccentricity (circular orbit)
i = deg2rad(105);   % Inclination (radians)
RAAN_0 = deg2rad(0); % Right Ascension of Ascending Node
omega = 0;          % Argument of periapsis
theta0 = linspace(0,2*pi,10);         % True anomaly
theta0(end) = [];

period = 2*pi*sqrt(a^3/mu_earth);

n = 1/(sqrt(a^3/mu_earth));
% [r, v] = kepler2cart(a, e, i, RAAN_0, omega, theta0);
% h = cross(r, v);  % Angular momentum vector

% Satellite position [Rsc attitude and control]
sat_x = a * (cos(RAAN_0) * cos(theta0) - sin(RAAN_0) * sin(theta0) * cos(i));
sat_y = a * (sin(RAAN_0) * cos(theta0) + cos(RAAN_0) * sin(theta0) * cos(i));
sat_z = a * sin(theta0) * sin(i);

sat_pos = [sat_x; sat_y; sat_z];


% Compute full orbit
theta = linspace(0, 2*pi, 100);
orbit_x = a * (cos(RAAN_0) * cos(theta) - sin(RAAN_0) * sin(theta) * cos(i));
orbit_y = a * (sin(RAAN_0) * cos(theta) + cos(RAAN_0) * sin(theta) * cos(i));
orbit_z = a * sin(theta) * sin(i);

%% SENSOR CONE
direction = sat_pos / norm(sat_pos);  % Normalized sensor direction
cone_half = deg2rad(10);  % Half-cone angle
cone_height = 50;      % Cone height ¡
% 
% psi = linspace(0, 2*pi, 100);
% H_cone = linspace(0,cone_height,100);
% 
% [H_grid, Psi_grid] = meshgrid(H_cone,psi);
% 
% % Cartesian coordinates of a point on the cone surface in LOCAL FRAME
% x_local = H_grid.*sin(cone_half).*cos(Psi_grid);
% y_local = H_grid.*sin(cone_half).*sin(Psi_grid);
% z_local = H_grid.*cos(cone_half);
% 
% xl_cone = x_local(:);
% yl_cone = y_local(:);
% zl_cone = z_local(:);
% 
% cone_local = [xl_cone, yl_cone, zl_cone]';
% 
% global cone_meshgrid;
% cone_meshgrid.X = x_local;
% cone_meshgrid.Y = y_local;
% cone_meshgrid.Z = z_local;
% cone_meshgrid.points = cone_local;
% 
% % Computation of the rotation matrix !!!!
% u = cross([0; 0; 1],direction);
% u = u/norm(u);
% v = cross(direction,u);
% v = v/norm(v);
% 
% R = [u, v, direction];
% cone_points = R*cone_meshgrid.points + sat_pos;
% 
% % Reshape the points again into a grid
% xc_cone = reshape(cone_points(1,:),100,100); % the 100 represent the number of thetas given and the number of rs given
% yc_cone = reshape(cone_points(2,:),100,100);
% zc_cone = reshape(cone_points(3,:),100,100);

%% CENTER OF GRID
[incli,azim,rad]=meshgrid(linspace(inclination_min_rad, inclination_max_rad, 20),linspace(azimuth_min_rad, azimuth_max_rad, 20),linspace(radius_min,radius_min,1));
xg = rad .* cos(incli) .* cos(azim);
yg = rad .* cos(incli) .* sin(azim);
zg = rad .* sin(incli);

incli = reshape(incli,[],1);
azim = reshape(azim,[],1);
rad = reshape(rad,[],1);

xg = reshape(xg,[],1);
yg = reshape(yg,[],1);
zg = reshape(zg,[],1);

%% Simulation
J2 = 1.08263e-3;
RAAN_dot = (-3*pi*J2*Re^2*cos(i)/a^2)/period; % [rad/s]
Tsim = 2*pi/RAAN_dot/2; %[s]
%Tsim = 1000; 
t = 0;
%deltaT = (deg2rad(1.8)/(pi/period))/5; % [s]
deltaT = period/(2*pi)*0.00025*0.5; % COMPUTED IN THE IPAD!!! SAME AS NUMSATS
max_dist = 50;

figure;
hold on
[x,y,z]=sphere(100);
x=x*Re;
y=y*Re;
z=z*Re;
surf(x,y,z,"EdgeColor","none","FaceAlpha",.2,"FaceColor",[0.3010 0.7450 0.9330]);

s1 = scatter3(xg,yg,zg,'b','Marker','o','MarkerFaceColor','flat');
hold on;
s2 = scatter3([],[],[],'g','Marker','o','MarkerFaceColor','flat');
s3 = scatter3(sat_x,sat_y,sat_z,'r','Marker','o','MarkerFaceColor','flat');
axis equal;

%ptraza=plot3([],[],[],'Color','k');
%traza=zeros(3,500);
axis equal;
view(3);


visited = zeros(size(xg));
P =0;
b = 0;
while t<Tsim 
    am = tic;
    for k = 1:length(theta0)
        theta_in = theta0(k);
        [sat_x, sat_y, sat_z] = updateSatellite(t,i,a,theta_in,RAAN_dot);
        index = abs(sat_x-xg)<max_dist & abs(sat_y-yg)<max_dist & abs(sat_z-zg)<max_dist;
        possiblegrid=[xg(index) yg(index) zg(index)];
        for j = 1:size(possiblegrid,1)
                x = possiblegrid(j,1)-sat_x; y =possiblegrid(j,2)-sat_y; z=possiblegrid(j,3)-sat_z;
                d = [x y z];
                d_rel = d/norm(d);
                angle = atan2(norm(cross(d_rel,[sat_x, sat_y, sat_z])),dot(d_rel,[sat_x, sat_y, sat_z]));
                if angle<=cone_half && angle>=-cone_half && norm(d) < max_dist
                    idx2= xg==possiblegrid(j,1) & yg ==possiblegrid(j,2) & zg ==possiblegrid(j,3);
                    visited(idx2) = visited(idx2) + 1;
                end
        end
        sats_x(k) = sat_x; sats_y(k) = sat_y; sats_z(k) = sat_z;
    end
    if mod(P,10000) == 0
       set(s2, 'XData',xg(visited>0), 'YData',yg(visited>0), 'ZData',zg(visited>0));
       set(s3, 'XData',sats_x, 'YData',sats_y, 'ZData',sats_z)
       %delete(ptraza)
       %traza(:,1:end-1)= traza(:,2:end);
       %traza(:,end)=[sats_x(end); sats_y(end); sats_z(end)];
       %ptraza=plot3(traza(1,:),traza(2,:),traza(3,:),'color','k');
        drawnow;
        frame = getframe(gcf);
        writeVideo(videoFile, frame);
        disp(b)

    disp(100*t/Tsim)
    end
    t = t+deltaT;
    b = toc(am);
    P = P+1;
end

%% Plot


