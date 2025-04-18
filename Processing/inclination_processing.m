%% Inclination processing with two RAANS
clear all; close all; clc;

% RAAN0 = 310 deg
n_files = 5;
percentage = [];
percentage_total = [];

for i = 1:n_files
    % Load file dynamically
    load(sprintf('full_inclination%d.mat', i));
    
    % Compute visited and revisited
    visited = sum(visited > 0, 2);
    revisited = sum(revisited > 0, 2);
    
    % Compute percentages
    perc = revisited ./ visited;
    perc_total = revisited ./ n_debris;
    
    % Store results
    percentage = [percentage, perc];
    percentage_total = [percentage_total, perc_total];
end

% Convert to percentage
percentage = 100 * percentage;
percentage_total = 100 * percentage_total;

% Compute mean and std
percentage_mean = mean(percentage, 2);
percentage_total_mean = mean(percentage_total, 2);

percentage_std = std(percentage, 0, 2);
percentage_total_std = std(percentage_total, 0, 2);

% RAAN0 = 0 deg
n_files = 5;
percentage0 = [];
percentage_total0 = [];

for i = 1:n_files
    % Load file dynamically
    load(sprintf('full0_inclination%d.mat', i));
    
    % Compute visited and revisited
    visited = sum(visited > 0, 2);
    revisited = sum(revisited > 0, 2);
    
    % Compute percentages
    perc = revisited ./ visited;
    perc_total = revisited ./ n_debris;
    
    % Store results
    percentage0 = [percentage0, perc];
    percentage_total0 = [percentage_total0, perc_total];
end

% Convert to percentage
percentage0 = 100 * percentage0;
percentage_total0 = 100 * percentage_total0;

% Compute mean and std
percentage_mean0 = mean(percentage0, 2);
percentage_total_mean0 = mean(percentage_total0, 2);

percentage_std0 = std(percentage0, 0, 2);
percentage_total_std0 = std(percentage_total0, 0, 2);

% Figure 1 % of visited debris
figure;
errorbar(i_variation,percentage_mean,percentage_std,'o', ...                  % circle marker
         'MarkerSize', 7,'LineWidth',1.5);
hold on;
errorbar(i_variation,percentage_mean0,percentage_std0,'o', ...                  % circle marker
         'MarkerSize', 7,'LineWidth',1.5);
set(gca, 'FontSize',12,'FontName','Arial');
xlabel('Inclination [°]','FontSize',12,'FontName','Arial');
ylabel('Percentage of visited debris revisited [%]','FontSize',12,'FontName','Arial');
legend('\Omega_0 = 315° ','\Omega_0 = 0°')

% Figure 2 % of total debris
figure;
errorbar(i_variation,percentage_total_mean,percentage_total_std,'o', ...                  % circle marker
         'MarkerSize', 7,'LineWidth',1.5);
hold on;
errorbar(i_variation,percentage_total_mean0,percentage_total_std0,'o', ...                  % circle marker
         'MarkerSize', 7,'LineWidth',1.5);
set(gca, 'FontSize',12,'FontName','Arial');
xlabel('Inclination [°]','FontSize',12,'FontName','Arial');
ylabel('Percentage of total debris revisited [%]','FontSize',12,'FontName','Arial');
legend('\Omega_0 = 315° ','\Omega_0 = 0°','Location','best')


