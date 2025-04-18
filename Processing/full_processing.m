%% Revisit time processing
clear all; close all; clc;
% #1
load("full_simulation11_20min.mat")
%perc = sum(revisited>0)/sum(visited>0)*100;

temps = zeros(1,n_debris);

for u = 1:n_debris
    t_deb = TimeObserved{u};
    if ~isempty(t_deb)
        temps(u) = min(t_deb);
    end
end
temps(temps(:)==0) = [];

temps = sort(temps);

maxtime = max(temps);

timeInstants = 0:maxtime;
cumulativePercentage = zeros(size(timeInstants));

for idx = 1:length(timeInstants)
    currentTime = timeInstants(idx);
    
    % Count objects detected (nonzero) on or before the current time
    detectedCount = sum(temps > 0 & temps <= currentTime);
    cumulativePercentage(idx) = (detectedCount / n_debris) * 100;
end

% #2
load("full_simulation12_20min.mat")
temps2 = zeros(1,n_debris);
for u = 1:n_debris
    t_deb = TimeObserved{u};
    if ~isempty(t_deb)
        temps2(u) = min(t_deb);
    end
end
temps2(temps2(:)==0) = [];
temps2 = sort(temps2);

maxtime2 = max(temps2);

timeInstants2 = 0:maxtime2;
cumulativePercentage2 = zeros(size(timeInstants2));

for idx = 1:length(timeInstants2)
    currentTime = timeInstants2(idx);
    % Count objects detected (nonzero) on or before the current time
    detectedCount2 = sum(temps2 > 0 & temps2 <= currentTime);
    cumulativePercentage2(idx) = (detectedCount2 / n_debris) * 100;
end

% #3
load("full_simulation13_20min.mat")
temps3 = zeros(1,n_debris);
for u = 1:n_debris
    t_deb = TimeObserved{u};
    if ~isempty(t_deb)
        temps3(u) = min(t_deb);
    end
end
temps3(temps3(:)==0) = [];
temps3 = sort(temps3);

maxtime3 = max(temps3);

timeInstants3 = 0:maxtime3;
cumulativePercentage3 = zeros(size(timeInstants3));

for idx = 1:length(timeInstants3)
    currentTime = timeInstants3(idx);
    % Count objects detected (nonzero) on or before the current time
    detectedCount3 = sum(temps3 > 0 & temps3 <= currentTime);
    cumulativePercentage3(idx) = (detectedCount3 / n_debris) * 100;
end

% #4
load("full_simulation14_20min.mat")
temps4 = zeros(1,n_debris);
for u = 1:n_debris
    t_deb = TimeObserved{u};
    if ~isempty(t_deb)
        temps4(u) = min(t_deb);
    end
end
temps4(temps4(:)==0) = [];
temps4 = sort(temps4);

maxtime4 = max(temps4);

timeInstants4 = 0:maxtime4;
cumulativePercentage4 = zeros(size(timeInstants4));

for idx = 1:length(timeInstants4)
    currentTime = timeInstants4(idx);
    % Count objects detected (nonzero) on or before the current time
    detectedCount4 = sum(temps4 > 0 & temps4 <= currentTime);
    cumulativePercentage4(idx) = (detectedCount4 / n_debris) * 100;
end

% #5
load("full_simulation15.mat")
temps5 = zeros(1,n_debris);
for u = 1:n_debris
    t_deb = TimeObserved{u};
    if ~isempty(t_deb)
        temps5(u) = min(t_deb);
    end
end
temps5(temps5(:)==0) = [];
temps5 = sort(temps5);

maxtime5 = max(temps5);

timeInstants5 = 0:maxtime5;
cumulativePercentage5 = zeros(size(timeInstants5));

for idx = 1:length(timeInstants5)
    currentTime = timeInstants5(idx);
    % Count objects detected (nonzero) on or before the current time
    detectedCount5 = sum(temps5 > 0 & temps5 <= currentTime);
    cumulativePercentage5(idx) = (detectedCount5 / n_debris) * 100;
end

% #6
load("full_simulation16_20min.mat")
temps6 = zeros(1,n_debris);
for u = 1:n_debris
    t_deb = TimeObserved{u};
    if ~isempty(t_deb)
        temps6(u) = min(t_deb);
    end
end
temps6(temps6(:)==0) = [];
temps6 = sort(temps6);

maxtime6 = max(temps6);

timeInstants6 = 0:maxtime6;
cumulativePercentage6 = zeros(size(timeInstants6));

for idx = 1:length(timeInstants6)
    currentTime = timeInstants6(idx);
    % Count objects detected (nonzero) on or before the current time
    detectedCount6 = sum(temps6 > 0 & temps6 <= currentTime);
    cumulativePercentage6(idx) = (detectedCount6 / n_debris) * 100;
end


%% OBTAINED VALUES REVISITED/VISITED FOR THE 20 MIN SIMULATIONS
load("percentages_full.mat")
clearvars -except timeInstants timeInstants2 timeInstants3 timeInstants4 timeInstants5 timeInstants6 cumulativePercentage cumulativePercentage2 cumulativePercentage3 cumulativePercentage4 cumulativePercentage5 cumulativePercentage6
timeVecs = {timeInstants, timeInstants2, timeInstants3, timeInstants4, timeInstants5, timeInstants6};
percVecs = {cumulativePercentage, cumulativePercentage2, cumulativePercentage3, cumulativePercentage4, cumulativePercentage5, cumulativePercentage6};

tMax = max(cellfun(@(v) max(v(:)), timeVecs));  

tMin = 0;  % Assuming you start at time 0
 
% Create a common time grid covering the full range from all datasets.

commonTime = tMin:1:tMax;  % Adjust the time step as needed
 
numDatasets = numel(timeVecs);

interpPerc = zeros(numDatasets, length(commonTime));
 
% Loop over each dataset: sort the vectors, use nearest value, and extend as needed

for i = 1:numDatasets

    % Sort the time vector and corresponding percentage vector

    [xiSorted, sortIdx] = sort(timeVecs{i});

    yiSorted = percVecs{i}(sortIdx);

    % Use the 'nearest' method to pick the closest value from the dataset

    % For points outside the range, we initially assign NaN.

    interpolated = interp1(xiSorted, yiSorted, commonTime, 'nearest', NaN);

    % For time instants beyond the maximum time in the dataset, 

    % set the value to the last recorded percentage.

    beyondIdx = commonTime > xiSorted(end);

    interpolated(beyondIdx) = yiSorted(end);

    % Optionally, for time instants before the first time point, 

    % assign the first value (if needed).

    beforeIdx = commonTime < xiSorted(1);

    interpolated(beforeIdx) = yiSorted(1);

    interpPerc(i,:) = interpolated;

end
 
% Compute the mean and standard deviation across all datasets for each common time instant
meanPercentage = mean(interpPerc);
stdPercentage  = std(interpPerc);  % Standard deviation computed along each time instant

commonTime = commonTime./(24*3600);

% Plot the mean cumulative percentage with error bars representing one standard deviation
figure;
errorbar(commonTime/(3600/24), meanPercentage, stdPercentage, 'LineWidth', 2);
plot(commonTime, meanPercentage, 'LineWidth', 2);
set(gca, 'FontSize', 12, 'FontName', 'Arial');
xlabel('Time [days]','FontSize', 12, 'FontName', 'Arial');
ylabel('Percentage of debris detected [%]','FontSize', 12, 'FontName', 'Arial');
%title('Debris detected by RADAR constellation','FontSize', 16, 'FontName', 'Arial');
%xlim([0,400]);
grid on;


x = commonTime(:);
y = meanPercentage(:);

% NON-LINEAR
modelo = @(beta,x) beta(1)*log(1+x./beta(2));
beta0 = [1; 1];
beta_est = nlinfit(x,y,modelo,beta0);

x_fit = linspace(0,max(x)*100,500);
y_fit = beta_est(1)*log(1+x_fit./beta_est(2));

figure;
plot(x,y,'bo','DisplayName','Data');
hold on;
plot(x_fit,y_fit,'Color',[0 0 0],'LineStyle','--','DisplayName','Fitted');

ylim([0 80]);
xlim([0 10000])

xline(5*365,':','Color',[0 0 0]);
ylims = ylim;  % get current y-axis limits
text(5*365 + 150, ylims(2)*0.80, '5 years', 'Rotation', 90, ...
    'HorizontalAlignment', 'left', 'FontSize', 12);


xlabel('Time [days]','FontSize', 12, 'FontName', 'Arial');
ylabel('Percentage of debris catalogued [%]','FontSize', 12, 'FontName', 'Arial');

legend('Data','Fitted','Location','southeast');
%%
log_x = log(x+eps);

coeffs = polyfit(log_x,y,1);
a = coeffs(1);
b = coeffs(2);

x_extrap = linspace(max(x),max(x)*10,1000);

y_extrap = a*log(x_extrap)+b;

figure;
plot(x,y,'bo');
hold on;

x_fit = linspace(min(x),max(x),1000);
y_fit = a*log(x_fit)+b;
plot(x_fit,y_fit,'r-','LineWidth',2);
hold on;

plot(x_extrap,y_extrap, 'k--','LineWidth',2)
hold on;

legend('Original','Adjust','Extrapolated')

xlim([min(x) max(x_extrap)])



%% ASSUMING LINEAR

[fitObj,gof] = fit(commonTime', meanPercentage', 'poly1');  % Can try 'poly1', 'exp1', etc.
disp(gof.rsquare);  % This is R²

% Find when 80% is covered
target_y = 80;

a = fitObj.p1;  % slope
b = fitObj.p2;  % intercept

x_vals = linspace(0, 50, 500);  % Change range based on your use case

% Compute y
y_vals = a * exp(b * x_vals);

% Plot it
figure;
plot(x_vals, y_vals, 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('y');
title('Exponential Curve: y = a \cdot e^{b \cdot x}');
grid on;


a = fitObj.p1;  % slope
b = fitObj.p2;  % intercept

target_x = (target_y - b) / a;

fprintf('x for y = %.2f is approximately %.4f\n years', target_y, target_x/365);

% Find how much is covered in 5 years
x_target = 5*365;

y_target = x_target*a + b;

fprintf('x for y = %.2f is approximately %.4f\n years', y_target, x_target/365);


%%
% percentage = [17.2308 13.4286 13.1195 16.0526 17.1512 15.3614];
% 
% mean_perc = mean(percentage);
% 
% std_perc = std(percentage);
