%% RADAR PROCESSING
clear all; close all; clc;

load("radar_full1.mat")
 
% Find the maximum detection time (ignoring zeros)
maxTime = max(time);
 
% Create a vector of time instants (adjust the step if needed)
timeInstants = 0:maxTime;  
cumulativePercentage = zeros(size(timeInstants));
 
% Compute the cumulative percentage for each time instant
for idx = 1:length(timeInstants)
    currentTime = timeInstants(idx);
    % Count objects detected (nonzero) on or before the current time
    detectedCount = sum(time > 0 & time <= currentTime);
    cumulativePercentage(idx) = (detectedCount / n_debris) * 100;
end
%%
load("radar_full2.mat")
 
% Find the maximum detection time (ignoring zeros)
maxTime2 = max(time);
 
% Create a vector of time instants (adjust the step if needed)
timeInstants2 = 0:maxTime2;  
cumulativePercentage2 = zeros(size(timeInstants2));
 
% Compute the cumulative percentage for each time instant
for idx = 1:length(timeInstants2)
    currentTime2 = timeInstants2(idx);
    % Count objects detected (nonzero) on or before the current time
    detectedCount2 = sum(time > 0 & time <= currentTime2);
    cumulativePercentage2(idx) = (detectedCount2 / n_debris) * 100;
end

load("radar_full3.mat")
 
% Find the maximum detection time (ignoring zeros)
maxTime3 = max(time);
 
% Create a vector of time instants (adjust the step if needed)
timeInstants3 = 0:maxTime3;  
cumulativePercentage3 = zeros(size(timeInstants3));
 
% Compute the cumulative percentage for each time instant
for idx = 1:length(timeInstants3)
    currentTime3 = timeInstants3(idx);
    % Count objects detected (nonzero) on or before the current time
    detectedCount3 = sum(time > 0 & time <= currentTime3);
    cumulativePercentage3(idx) = (detectedCount3 / n_debris) * 100;
end

load("radar_full4.mat")
 
% Find the maximum detection time (ignoring zeros)
maxTime4 = max(time);
 
% Create a vector of time instants (adjust the step if needed)
timeInstants4 = 0:maxTime4;  
cumulativePercentage4 = zeros(size(timeInstants4));
 
% Compute the cumulative percentage for each time instant
for idx = 1:length(timeInstants4)
    currentTime4 = timeInstants4(idx);
    % Count objects detected (nonzero) on or before the current time
    detectedCount4 = sum(time > 0 & time <= currentTime4);
    cumulativePercentage4(idx) = (detectedCount4 / n_debris) * 100;
end

load("radar_full5.mat")
 
% Find the maximum detection time (ignoring zeros)
maxTime5 = max(time);
 
% Create a vector of time instants (adjust the step if needed)
timeInstants5 = 0:maxTime5;  
cumulativePercentage5 = zeros(size(timeInstants5));
 
% Compute the cumulative percentage for each time instant
for idx = 1:length(timeInstants5)
    currentTime5 = timeInstants5(idx);
    % Count objects detected (nonzero) on or before the current time
    detectedCount5 = sum(time > 0 & time <= currentTime5);
    cumulativePercentage5(idx) = (detectedCount5 / n_debris) * 100;
end

%%
load("percentages_radar.mat")
clearvars -except timeInstants timeInstants2 timeInstants3 timeInstants4 timeInstants5 cumulativePercentage cumulativePercentage2 cumulativePercentage3 cumulativePercentage4 cumulativePercentage5
timeVecs = {timeInstants, timeInstants2, timeInstants3, timeInstants4, timeInstants5};
percVecs = {cumulativePercentage, cumulativePercentage2, cumulativePercentage3, cumulativePercentage4, cumulativePercentage5};

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
y = meanPercentage(:)./100;

% Make sure x > -1 for log(x + 1)
validIdx = x > -1;
x = x(validIdx);
y = y(validIdx);

% Define the model: y = a*log(x + 1)
ft = fittype('a*log(x + 1)', 'independent', 'x', 'coefficients', 'a');
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.StartPoint = 1;

[fitObj, gof] = fit(x, y, ft, opts);

