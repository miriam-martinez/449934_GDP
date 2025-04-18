clear all; close all; clc;
% Constants and Earth parameters
Re = 6378;                % Earth's radius in km
mu_earth = 398600;        % Earth's gravitational parameter [km^3/s^2]
J2 = 1.08263e-3; 

%% Generate Sensor (CAMERA)
cone_half = deg2rad(6);  % Half-cone angle
range = 800-600;
maxrange = 1000;
minrange = 100;
%x_cone = minrange*tan(cone_half);

%% Generate Constellation
% Compute Last Layer
H_sat = [657.96 742.04];
%H_sat = 657.96;
%H_sat = 742.04; % [km]

% Satellites definition
a = H_sat + Re;
period = 2*pi*sqrt(a.^3/mu_earth);
inc = [deg2rad(98.015) deg2rad(98.357)];
%inc = deg2rad(98.015);
n = 1./(sqrt(a.^3/mu_earth));

% Compute number of satellites
%RAAN_dot = -3*pi*J2*Re^2*cos(inc)/(H_sat+Re)^2; % [rad/s]
RAAN_dot = 1.991e-7; % [rad/s]
%tscan = 2*pi/abs(RAAN_dot)*period;
t_eclipse = 6.64; % [min]
theta_ecl = 2*pi./period.*t_eclipse*60; % [rad]


%% Generate Debris
n_debris = 100;

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

%% SATS
n_sats = 27*2;
deltaTheta=deg2rad(0.685);
for iterative=1:2:n_sats
    theta0_orbit(iterative)=(iterative-1)*2*pi/(n_sats);
    theta0_orbit(iterative+1)=theta0_orbit(iterative)+deltaTheta;
end
theta0 = [theta0_orbit, theta0_orbit];
n_sats = size(H_sat,2)*n_sats;
coord_sat = zeros(3,n_sats);

axcone=zeros(2,1);


%% Simulation
Tsim = 7*24*60*60;
t = 0;
deltaT = 1;

visited = zeros(1,n_debris);
timing=cell(1,n_debris);
P = 0;
b = 0;
LastTimeVisitedDebris=-inf(n_debris,1);
FirstTimeVisitedDebris=-inf(n_debris,1);
while t<Tsim 
    am = tic;
     [coord_sat,coord_deb,axcone,eclipse]=updateOrbitalElement(a, RAAN_dot, t, inc, n, theta0, RAAN0d, RAANd_dot, a_d, i_d,  theta0d, n_d, theta_ecl);
    for k = 1:1:n_sats
        if eclipse(k)
            continue;
        end
        index = abs(coord_sat(1,k)-coord_deb(1,:))<maxrange & abs(coord_sat(2,k)-coord_deb(2,:))<maxrange & abs(coord_sat(3,k)-coord_deb(3,:))<maxrange;
        possibledebris = coord_deb(:,index);
        for j = 1:size(possibledebris,2)
                d = possibledebris(:,j) - coord_sat(:,k);
                angle = atan2(norm(cross(d,axcone)),dot(d,axcone));
                if angle<=cone_half && angle>=-cone_half  && norm(d)<maxrange && norm(d)>minrange
                    idx2 = all(coord_deb(:,:) == possibledebris(:,j));
                    
                    if t-LastTimeVisitedDebris(idx2)<1.1*deltaT
                            LastTimeVisitedDebris(idx2)=t;
                    else
                        LastTimeVisitedDebris(idx2)=t;
                        FirstTimeVisitedDebris(idx2)=t;
                    end

                    if t-FirstTimeVisitedDebris(idx2)==10
                        visited(idx2) = visited(idx2) + 1;
                        timing{idx2}=[t,timing{idx2}];
                        continue;
                    end
                end
        end
    end
    if mod(P,5000) == 0
        disp(sum(visited>0)+" of "+n_debris)
        disp("t="+t/60/60+"h")
    end
    t = t+deltaT;
    b = toc(am);
    P = P+1;
end



function [coord_sat,coord_deb,axcone,eclipse]=updateOrbitalElement( a, RAAN_dot, t, inc, n, theta0, RAAN0d, RAANd_dot, a_d, i_d, theta0d, n_d,theta_ecl)
%global coord_sat a RAAN_dot t inc n theta0 RAAN0d RAANd_dot a_d i_d coord_deb theta0d n_d
RAAN = deg2rad(0) + RAAN_dot*t;
cosRAAN = cos(RAAN)*ones(1,length(theta0)/2);
sinRAAN = sin(RAAN)*ones(1,length(theta0)/2);

theta = reshape(theta0,54,2) + n*t;
costheta = cos(theta);
sintheta = sin(theta);

eclipse1 = zeros(1,length(theta));
eclipse2 = zeros(1,length(theta));

eclipse1(mod(theta(:,1)',2*pi)>0 & mod(theta(:,1)',2*pi)<theta_ecl(1)) = 1;
eclipse2(mod(theta(:,2)',2*pi)>0 & mod(theta(:,2)',2*pi)<theta_ecl(2)) = 1;
eclipse = [eclipse1 eclipse2];

sini = sin(inc);
cosi = cos(inc);

coord_sat(1,1:54) = a(1).* (cosRAAN.* costheta(:,1)' - sinRAAN.* sintheta(:,1)'.* cosi(1));
coord_sat(2,1:54) = a(1).* (sinRAAN.* costheta(:,1)' + cosRAAN.* sintheta(:,1)'.* cosi(1));
coord_sat(3,1:54) = a(1).* sintheta(:,1)'.* sini(1);

coord_sat(1,55:108) = a(2).* (cosRAAN.* costheta(:,2)' - sinRAAN.* sintheta(:,2)'.* cosi(2));
coord_sat(2,55:108) = a(2).* (sinRAAN.* costheta(:,2)' + cosRAAN.* sintheta(:,2)'.* cosi(2));
coord_sat(3,55:108) = a(2).* sintheta(:,2)'.* sini(2);

axcone=[sinRAAN(1),cosRAAN(1),0];

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