%% Full constellation
clear; close all; clc;
rng(0)
% Constants and Earth parameters
Re = 6378;                % Earth's radius in km
mu_earth = 398600;        % Earth's gravitational parameter [km^3/s^2]
J2 = 1.08263e-3;

%% Generate Sensor (RADAR)
cone_half_rad = deg2rad(10);  % Half-cone angle
range_rad = 800-600;
maxrange_rad = 50;
minrange_rad = 5;
layers_rad = ceil(200/(maxrange_rad-minrange_rad));
min_r_rad = maxrange_rad-range_rad/layers_rad; 
x_cone_rad = min_r_rad*tan(cone_half_rad);

% Compute Layer
H_sat_rad = 600-min_r_rad + ([1:layers_rad]-1)*(maxrange_rad-min_r_rad);
hip_rad = Re + min_r_rad + H_sat_rad;

% Satellites definition
t_eclipse_rad = 30; % [min]
a_rad = H_sat_rad + Re;
period_rad = 2*pi*sqrt(a_rad.^3/mu_earth);
theta_ecl_rad = 2*pi./period_rad.*t_eclipse_rad*60;

%% Generate Sensor (CAMERA)
cone_half_cam = deg2rad(6);
maxrange_cam = 1000;
minrange_cam = 100;

%% Generate Debris
n_debris = 1000;

low_h = Re + 600;
up_h = Re + 800;
a_d = low_h+(up_h-low_h)*rand(1,n_debris);
n_d = 1./(sqrt(a_d.^3/mu_earth));

low_i = 60;
up_i =105;
i_d =deg2rad(low_i+(up_i-low_i)*rand(1,n_debris));

theta0d = (2*pi)*rand(1,n_debris);

RAAN0d = (2*pi)*rand(1,n_debris);
RAANd_dot = -3*pi*J2*Re^2*cos(i_d)./a_d.^2; % [rad/s]

coord_deb = zeros(3,n_debris);

%% Generate Constellation (RADAR)
inc_rad = deg2rad(80);
RAAN_0_rad = 0;
raan_dif_rad = 2*asin(x_cone_rad./hip_rad); % [rad]
n_rad = 1./(sqrt(a_rad.^3/mu_earth));
RAAN_dot_rad = -3*pi*J2*Re^2*cos(inc_rad)./(H_sat_rad+Re).^2; % [rad/s]
RAAN0_rad = 0;

orbit_rad = raan_dif_rad./abs(RAAN_dot_rad); % [orbit]
delta_theta_rad = 2*pi*orbit_rad; % [rad] 
n_sats_rad = ceil(2*pi./delta_theta_rad);
RAAN_dot_rad_concat=[];
theta0_rad_concat=[];
for idx11=1:layers_rad
    RAAN_dot_rad_concat=[RAAN_dot_rad_concat,repmat(RAAN_dot_rad(idx11),1,n_sats_rad(idx11))];
    theta0_rad = linspace(0,2*pi,n_sats_rad(idx11));
    theta0_rad_concat=[theta0_rad_concat, theta0_rad];
    theta0_rad(end) = [];
end
coord_sat_rad = zeros(3,sum(n_sats_rad));

%% Generate Constellation (OPTICAL)
H_sat_cam = [657.96 742.04];
a_cam = H_sat_cam + Re;
period_cam = 2*pi*sqrt(a_cam.^3/mu_earth);

inc_cam = [deg2rad(98.015) deg2rad(98.357)];
n_cam = 1./(sqrt(a_cam.^3/mu_earth));
RAAN_dot_cam = 1.991e-7;
RAAN0_cam = 0;

t_eclipse_cam = 6.64; % [min]
theta_ecl_cam = 2*pi./period_cam.*t_eclipse_cam*60; % [rad]

n_sats_cam = 27*2;
deltaTheta=deg2rad(0.685);
for iterative=1:2:n_sats_cam
    theta0_orbit(iterative)=(iterative-1)*2*pi/(n_sats_cam);
    theta0_orbit(iterative+1)=theta0_orbit(iterative)+deltaTheta;
end
theta0_cam = [theta0_orbit, theta0_orbit];
n_sats_cam = size(H_sat_cam,2)*n_sats_cam;
coord_sat_cam = zeros(3,n_sats_cam);

axcone=zeros(2,1);


%% Simulation
%Tsim = 3.5*30*24*3600; % [s]
Tsim = 1*24*3600; % [s]
t = 0;
%deltaT = period_rad(1)/(2*pi)*0.00025*0.5; % COMPUTED IN THE IPAD!! SAME AS NUMSATS
deltaT = 0.125;

% figure;
% s1 = scatter3([],[],[],'b','Marker','.');
% hold on;
% s2 = scatter3([],[],[],'g','Marker','o','MarkerFaceColor','flat');
% s3 = scatter3([],[],[],'r','Marker','o','MarkerFaceColor','flat');
% axis equal;


visited = zeros(1,n_debris);
revisited = zeros(1,n_debris);
timing=cell(1,n_debris);
TimeObserved = cell(1,n_debris);
P = 0;
b = 0;
LastTimeVisitedDebris=-inf(n_debris,1);
FirstTimeVisitedDebris=-inf(n_debris,1);
lasttime = -inf(n_debris,1);
firsttime = -inf(n_debris,1);
while t<Tsim 
    am = tic;
     [coord_sat_rad,coord_deb,eclipse_rad]=updateOrbitalElement(a_rad, RAAN0_rad, RAAN_dot_rad_concat, t, inc_rad, n_rad, theta0_rad_concat, RAAN0d, RAANd_dot, a_d, i_d,  theta0d, n_d, theta_ecl_rad,1);
    for k = 1:1:sum(n_sats_rad)
        if eclipse_rad(k)
            continue;
        end
        index = abs(coord_sat_rad(1,k)-coord_deb(1,:))<maxrange_rad & abs(coord_sat_rad(2,k)-coord_deb(2,:))<maxrange_rad & abs(coord_sat_rad(3,k)-coord_deb(3,:))<maxrange_rad;
        possibledebris = coord_deb(:,index);
        for j = size(possibledebris,2):-1:1
                d = possibledebris(:,j) - coord_sat_rad(:,k);
                angle = atan2(norm(cross(d,coord_sat_rad(:,k))),dot(d,coord_sat_rad(:,k)));
                if angle<=cone_half_rad && angle>=-cone_half_rad  && norm(d)<=maxrange_rad && norm(d)>=min_r_rad
                    idx2 = all(coord_deb(:,:) == possibledebris(:,j));
                    
                    if t-LastTimeVisitedDebris(idx2)<=deltaT
                        LastTimeVisitedDebris(idx2)=t;
                    else
                        LastTimeVisitedDebris(idx2)=t;
                        FirstTimeVisitedDebris(idx2)=t;
                    end

                    if t-FirstTimeVisitedDebris(idx2)>=0.125
                        visited(idx2) = visited(idx2)+1;
                        %coord_deb_cam = [];
                        for time = t:1:(t+(1*20*60)) 
                            [coord_sat_cam,coord_deb_cam,eclipse_cam,axcone]=updateOrbitalElement(a_cam, RAAN0_cam, RAAN_dot_cam, time, inc_cam, n_cam, theta0_cam, RAAN0d(idx2), RAANd_dot(idx2), a_d(idx2), i_d(idx2),  theta0d(idx2), n_d(idx2), theta_ecl_cam,0);
                            for nsc = 1:1:n_sats_cam
                                if eclipse_cam(nsc)
                                    continue
                                end
                                d2 = coord_deb_cam - coord_sat_cam(:,nsc);
                                angle2 = atan2(norm(cross(d2,axcone)),dot(d2,axcone));
                                if angle2<=cone_half_cam && angle2>=-cone_half_cam && norm(d2)<=maxrange_cam && norm(d2)>=minrange_cam
                                    
                                    if time-lasttime(idx2) <1.1
                                        lasttime(idx2)  = time;
                                    else
                                        lasttime(idx2)  = time;
                                        firsttime(idx2)  = time;
                                    end

                                    if time-firsttime(idx2)==10
                                        revisited(idx2) = revisited(idx2) + 1;
                                        fprintf("he entrado \n");
                                        TimeObserved{idx2} = [time, TimeObserved{idx2}];
                                        break
                                    end
                                end
                            end
                        end

                    end
                end
        end
    end
    if mod(P,50000) == 0
       % set(s1, 'XData',coord_deb(1,visited==0), 'YData',coord_deb(2,visited==0), 'ZData',coord_deb(3,visited==0));
       % set(s2, 'XData',coord_deb(1,visited>0), 'YData',coord_deb(2,visited>0), 'ZData',coord_deb(3,visited>0));
       % set(s3, 'XData',coord_sat(1,:), 'YData',coord_sat(2,:), 'ZData',coord_sat(3,:))
       %  drawnow;
        %disp(100*t/Tsim+"%")
        disp("Debris remaining: " + length(LastTimeVisitedDebris))
        %disp("Time remaining: "+b*(Tsim-t)/deltaT/60)
    end
    t = t+deltaT;
    b = toc(am);
    P = P+1;
end

nonVisited=length(LastTimeVisitedDebris)/n_debris;
 


function [coord_sat,coord_deb,eclipse,axcone]=updateOrbitalElement(a,RAAN0, RAAN_dot, t, inc, n, theta0, RAAN0d, RAANd_dot, a_d, i_d, theta0d, n_d,theta_ecl,flag)
RAAN = RAAN0 + RAAN_dot.*t;
if flag == 1
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
else
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
end

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