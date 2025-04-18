%% Optical Number of Satellites
clear all; close all; clc;
% Why pairs and not triplets?
load("lowjoined_optical.mat");
x = 1:1:n_sats;
mean_percentage = mean(percentage*100);
std_percentage = std(percentage*100);
x = categorical(x);

bar(x,mean_percentage)
hold on;
errorbar(x,mean_percentage,std_percentage)

set(gca, 'FontSize', 14, 'FontName', 'Arial');
xlabel('# satellites','FontSize', 14, 'FontName', 'Arial');
ylabel('Pecentage of debris [%]','FontSize', 14, 'FontName', 'Arial')


%{
% PAIRS ANALYSIS
load("upheight_numsats_optical.mat");
x = 1:1:n_sats/2;
mean_percentage = mean(percentage*100);
std_percentage = std(percentage*100);
x = categorical(x);

bar(x,mean_percentage)
hold on;
errorbar(x,mean_percentage,std_percentage)

set(gca, 'FontSize', 14, 'FontName', 'Arial');
xlabel('Pairs of satellites','FontSize', 14, 'FontName', 'Arial');
ylabel('Pecentage of debris [%]','FontSize', 14, 'FontName', 'Arial')

%title('Debris detected by Plane I','FontSize', 16, 'FontName', 'Arial')
%}


