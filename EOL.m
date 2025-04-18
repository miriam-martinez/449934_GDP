%% Low orbit EoL

% Re = 6378; % [km]
% Mu = 398600.4418; 
% a = 400 + Re;
% e = 0;
% Cd = 2.2;
% m0 = 29.16; % [kg]
% A = (0.06+(2.273*0.239));
% Ap = 0;
% F107 = 30;
% 
% [P,t] = computeOrbitalDecay(a,e,A,Cd,m0,F107,Ap);
% 
% 
% t = t./(3600*24);
% y = ((P./(2.*pi)).^2.*Mu).^(1/3)-Re;
% 
% plot(t,y,'k','linewidth',2);

%% FROM OPERATIONAL ORBIT TO PARKING ORBIT

clear all; close all; clc;
Re = 6378; % [km]
mu = 3.986*10^14; % [m^3/s^2]
r_ope = (Re+750)*1000;
N = [3, 54]; % satellites per plane


% ELECTRIC PROPULSION
Isp = 2250 * 9.81; % [m/s]
T = 105*10^(-6); % [N]
mf = [29.2 26.9]; % [kg]

Cd = 2.2;
A = (0.06+(2.273*0.239)); % [m^2]
H0_pos = [0 100 150 200 250 300 350 400 450 500 550 600 650 700 750]*1000;
    
delta_pos = [8.4 5.9 25.5 37.5 44.8 50.3 54.8 58.2 61.3 64.5 68.7 74.8 84.4 99.3 121];
rho_pos = [1.225 5.25e-7 1.73e-9 2.41e-10 5.97e-11 1.87e-11 6.66e-12 2.62e-12 1.09e-12 4.76e-13 2.14e-13 9.89e-14 4.73e-14 2.36e-14 1.24e-14];


for sat = 1:length(N)
    k =0;
    for h = 100:50:750
        r_parking = (h+Re)*1000;
        deltaM = 2*pi/N(sat);

        V_parking = sqrt(mu/(r_parking));
        V_ope = sqrt(mu/r_ope);
        
        deltaV = abs(V_ope-V_parking); % [m/s]
        
        n_ope = sqrt(mu/r_ope^3);
        n_parking = sqrt(mu/r_parking^3);
        
        a_phasing = (r_ope+r_parking)/2;
        
        deltaT = deltaM/(n_parking+n_ope);
        time_waiting = N(sat)*deltaT/(3600*24); % [days]
        
        mi_optical = mf(sat)*exp(deltaV/Isp); 
        mprop_optical = mi_optical - mf(sat);
        deltaTep_optical = ((mprop_optical*Isp)/T)/(3600*24*365);
        TOTAL_ep = (deltaTep_optical);
        
        TIME = time_waiting + TOTAL_ep;
        
        if TIME >=2
            continue
        end
        
        Hf = 0;
        H0 = h*1000; % [km]
        deltaH = (Hf-H0);
        [~,i] = min(abs(H0-H0_pos));
        delta = delta_pos(i)*1000;
        beta = 1/delta;
        rho_0 = rho_pos(i);    
        m = mf(sat); % [kg]
        
        deltaT = m/(beta*rho_0*A*Cd*sqrt(mu*Re))*(1-exp(beta*deltaH));
        deltaT_days = deltaT/(60*60*24);
        k = k+1;
        total_time(k,sat) = deltaT_days + TIME*365;
        time_prop(k,sat) = TIME;
        heights(k,sat) = h;
    
    end
end

plot(heights(:,1),total_time(:,1)/365,'Color','b','DisplayName','RADAR constellation')
hold on;
plot(heights(:,2),total_time(:,2)/365,'Color','r','DisplayName','Optical constellation')
ylim([0, 10])
xlabel('Orbtial Altitude [km]', 'FontSize', 14, 'FontName', 'Arial','LineWidth',2)
ylabel('Time to dacay [years]','FontSize', 14, 'FontName', 'Arial','LineWidth',2)
legend("show",'FontSize', 9, 'FontName', 'Arial')