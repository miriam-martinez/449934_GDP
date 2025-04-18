%% RAAN VARIATION vs CATALOGUED DEBRIS GRAPH (1 day)
clear all; close all; clc;
load("raan_variation1.mat")
visited1 = sum(visited>0,2);
revisited1 = sum(revisited>0,2);
percentage1 = revisited1./visited1;
percentage12 = revisited1./n_debris;

load("raan_variation2.mat")
visited2 = sum(visited>0,2);
revisited2 = sum(revisited>0,2);
percentage2 = revisited2./visited2;
percentage22 = revisited2./n_debris;

load("raan_variation3.mat")
visited3 = sum(visited>0,2);
revisited3 = sum(revisited>0,2);
percentage3 = revisited3./visited3;
percentage32 = revisited3./n_debris;

load("raan_variation4.mat")
visited4 = sum(visited>0,2);
revisited4 = sum(revisited>0,2);
percentage4 = revisited4./visited4;
percentage42 = revisited4./n_debris;

load("raan_variation5.mat")
visited5 = sum(visited>0,2);
revisited5 = sum(revisited>0,2);
percentage5 = revisited5./visited5;
percentage52 = revisited5./n_debris;



percentage = 100*[percentage1 percentage2 percentage3 percentage4 percentage5];
percentage_total = 100*[percentage12 percentage22 percentage32 percentage42 percentage52];


percentage_mean = mean(percentage,2);
percentage_total_mean = mean(percentage_total,2);

percentage_std = std(percentage,0,2);
percentage_total_std = std(percentage_total,0,2);

figure;
%scatter(RAAN0_variation,percentage_mean)
errorbar(RAAN0_variation,percentage_mean,percentage_std,'o', ...                  % circle marker
         'MarkerSize', 7,'LineWidth',1.5);
set(gca, 'FontSize',12,'FontName','Arial');
xlabel('\Omega_0[°]','FontSize',12,'FontName','Arial');
ylabel('Percentage of visited debris revisited [%]','FontSize',12,'FontName','Arial');

figure;
%scatter(RAAN0_variation,percentage_total_mean)
errorbar(RAAN0_variation,percentage_total_mean,percentage_total_std,'o', ...                  % circle marker
         'MarkerSize', 7,'LineWidth',1.5);
set(gca, 'FontSize',12,'FontName','Arial');
xlabel('\Omega_0[°]','FontSize',12,'FontName','Arial');
ylabel('Percentage of total debris revisited [%]','FontSize',12,'FontName','Arial');
