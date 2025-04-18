clear all; close all; clc;
% Constants and Earth parameters
Re = 6378;                % Earth's radius in km
mu_earth = 398600;        % Earth's gravitational parameter [km^3/s^2]
J2 = 1.08263e-3; 

%% Generate Sensor (RADAR)
cone_half = deg2rad(10);  % Half-cone angle
range = 800-600;
maxrange = 50;
minrange = 5;
layers = ceil(200/(maxrange-minrange));
min_r = maxrange-range/layers; 
x_cone = min_r*tan(cone_half);


%% Generate Constellation
% Compute Last Layer
H_sat = 600-min_r + ([1:layers]-1)*(maxrange-min_r);
hip = Re + min_r + H_sat;

% Satellites definition
a = H_sat + Re;
period = 2*pi*sqrt(a.^3/mu_earth);
inc = deg2rad(80);
n = 1./(sqrt(a.^3/mu_earth));
RAAN_dot = -3*pi*J2*Re^2*cos(inc)./(H_sat+Re).^2; % [rad/s]
RAAN0 = deg2rad(0);

% Eclipse definition
t_eclipse = 30; % [min]
theta_ecl = 2*pi./period.*t_eclipse*60; % [rad]

% Number of satellites 
raan_dif = 2*asin(x_cone./hip); % [rad]
orbit = raan_dif./abs(RAAN_dot); % [orbit]
delta_theta = 2*pi.*orbit; % [rad] 
n_sats = ceil(2*pi./delta_theta);
theta0_orbit = linspace(0,2*pi,n_sats(1)+1);
theta0_orbit(end) = [];

theta0 = [theta0_orbit, theta0_orbit, theta0_orbit, theta0_orbit, theta0_orbit];
n_sats = sum(n_sats);
coord_sat = zeros(3,n_sats);

%% Generate Debris
n_debris = 1000;

low_h = 600 + Re;
up_h = 800 + Re;
a_d = low_h+(up_h-low_h)*rand(1,n_debris);
n_d = 1./(sqrt(a_d.^3/mu_earth));

low_i = 60;
up_i =105;
i_d =deg2rad(low_i+(up_i-low_i)*rand(1,n_debris));

theta0d = (2*pi)*rand(1,n_debris);

RAAN0d = (2*pi)*rand(1,n_debris);
RAANd_dot = -3*pi*J2*Re^2*cos(i_d)./a_d.^2; % [rad/s]

coord_deb = zeros(3,n_debris);


%% Simulation
Tsim = 3*365*60*60*24;
t = 0;
deltaT = 0.125;

visited = zeros(1,n_debris);
timing=cell(1,n_debris);
P = 0;
b = 0;
LastTimeVisitedDebris=-inf(n_debris,1);
FirstTimeVisitedDebris=-inf(n_debris,1);
percentage = 0;
time = zeros(1,n_debris);
while percentage<0.9 & t<Tsim 
    am = tic;
     [coord_sat,coord_deb,eclipse]=updateOrbitalElement(a,RAAN0, RAAN_dot, t, inc, n, theta0, RAAN0d, RAANd_dot, a_d, i_d, theta0d, n_d, theta_ecl);
    for k = 1:1:n_sats
        if eclipse(k)
            continue;
        end
        index = abs(coord_sat(1,k)-coord_deb(1,:))<maxrange & abs(coord_sat(2,k)-coord_deb(2,:))<maxrange & abs(coord_sat(3,k)-coord_deb(3,:))<maxrange;
        possibledebris = coord_deb(:,index);
        for j = size(possibledebris,2):-1:1
                d = possibledebris(:,j) - coord_sat(:,k);
                angle = atan2(norm(cross(d,coord_sat(:,k))),dot(d,coord_sat(:,k)));
                if angle<=cone_half && angle>=-cone_half  && norm(d)<maxrange && norm(d)>minrange
                    idx2 = all(coord_deb(:,:) == possibledebris(:,j));
                    
                    if t-LastTimeVisitedDebris(idx2)<1.1*deltaT
                            LastTimeVisitedDebris(idx2)=t;
                    else
                        LastTimeVisitedDebris(idx2)=t;
                        FirstTimeVisitedDebris(idx2)=t;
                    end

                    if t-FirstTimeVisitedDebris(idx2)>=0.125
                        visited(idx2) = visited(idx2) + 1;
                        timing{idx2}=[t,timing{idx2}];
                        percentage = sum(visited>0)/n_debris;
                        if time(idx2) == 0
                            time(idx2) = t;
                        end
                        continue;
                    end
                end
        end
    end
    if mod(P,500000) == 0
        disp(sum(visited>0)+" of "+n_debris)
        disp("t="+t/60/60/24+"days")
        disp(percentage*100 + "%")
    end
    t = t+deltaT;
    b = toc(am);
    P = P+1;
end



function [coord_sat,coord_deb,eclipse]=updateOrbitalElement(a,RAAN0, RAAN_dot, t, inc, n, theta0, RAAN0d, RAANd_dot, a_d, i_d, theta0d, n_d,theta_ecl)
RAAN = RAAN0 + RAAN_dot.*t;
cosRAAN = cos(RAAN);
sinRAAN = sin(RAAN);

theta = reshape(theta0,3,5) + n.*t;
costheta = cos(theta);
sintheta = sin(theta);

eclipse1 = zeros(1,length(theta));
eclipse2 = zeros(1,length(theta));
eclipse3 = zeros(1,length(theta));
eclipse4 = zeros(1,length(theta));
eclipse5 = zeros(1,length(theta));

eclipse1(mod(theta(:,1)',2*pi)>0 & mod(theta(:,1)',2*pi)<theta_ecl(1)) = 1;
eclipse2(mod(theta(:,2)',2*pi)>0 & mod(theta(:,2)',2*pi)<theta_ecl(2)) = 1;
eclipse3(mod(theta(:,2)',2*pi)>0 & mod(theta(:,2)',2*pi)<theta_ecl(2)) = 1;
eclipse4(mod(theta(:,2)',2*pi)>0 & mod(theta(:,2)',2*pi)<theta_ecl(2)) = 1;
eclipse5(mod(theta(:,2)',2*pi)>0 & mod(theta(:,2)',2*pi)<theta_ecl(2)) = 1;
eclipse = [eclipse1 eclipse2 eclipse3 eclipse4 eclipse5];

sini = sin(inc);
cosi = cos(inc);

coord_sat(1,1:3) = a(1).* (cosRAAN(1).* costheta(:,1)' - sinRAAN(1).* sintheta(:,1)'.* cosi);
coord_sat(2,1:3) = a(1).* (sinRAAN(1).* costheta(:,1)' + cosRAAN(1).* sintheta(:,1)'.* cosi);
coord_sat(3,1:3) = a(1).* sintheta(:,1)'.* sini;

coord_sat(1,4:6) = a(2).* (cosRAAN(2).* costheta(:,2)' - sinRAAN(2).* sintheta(:,2)'.* cosi);
coord_sat(2,4:6) = a(2).* (sinRAAN(2).* costheta(:,2)' + cosRAAN(2).* sintheta(:,2)'.* cosi);
coord_sat(3,4:6) = a(2).* sintheta(:,2)'.* sini;

coord_sat(1,7:9) = a(3).* (cosRAAN(3).* costheta(:,3)' - sinRAAN(3).* sintheta(:,3)'.* cosi);
coord_sat(2,7:9) = a(3).* (sinRAAN(3).* costheta(:,3)' + cosRAAN(3).* sintheta(:,3)'.* cosi);
coord_sat(3,7:9) = a(3).* sintheta(:,1)'.* sini;

coord_sat(1,10:12) = a(4).* (cosRAAN(4).* costheta(:,4)' - sinRAAN(4).* sintheta(:,4)'.* cosi);
coord_sat(2,10:12) = a(4).* (sinRAAN(4).* costheta(:,4)' + cosRAAN(4).* sintheta(:,4)'.* cosi);
coord_sat(3,10:12) = a(4).* sintheta(:,4)'.* sini;

coord_sat(1,13:15) = a(5).* (cosRAAN(5).* costheta(:,5)' - sinRAAN(5).* sintheta(:,5)'.* cosi);
coord_sat(2,13:15) = a(5).* (sinRAAN(5).* costheta(:,5)' + cosRAAN(5).* sintheta(:,5)'.* cosi);
coord_sat(3,13:15) = a(5).* sintheta(:,5)'.* sini;


% DEBRIS
RAANd = RAAN0d + RAANd_dot*t;
cosRAANd = cos(RAANd);
sinRAANd = sin(RAANd);

thetad = theta0d + n_d*t;
costhetad = cos(thetad);
sinthetad = sin(thetad);

sinid = sin(i_d);
cosid = cos(i_d);

coord_deb(1,:) = a_d.* (cosRAANd.* costhetad - sinRAANd.* sinthetad.* cosid);
coord_deb(2,:) = a_d.* (sinRAANd.* costhetad + cosRAANd.* sinthetad.* cosid);
coord_deb(3,:) = a_d.* sinthetad.* sinid;

end