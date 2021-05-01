% Main function for calculating the number of neurons in culture
% Requires dependencies:
%   manualCount.m
%   parseCoord.m
%   countNeurons.m
%   TotalCellMatterArea.m
%   WeightedDarkness.m
%   countNeurons.mcorrect

clearvars;
close all;
% Clears workspace.
imname = 'richnhighres10.jpg';
% Name of image to be analysed.
I = imread(imname);
% Reads in image.
windowSize = 300;
% Size of the cropped image.
[ROI, xl, yl] = manualCount(imname, windowSize); % Manual counting
stats = countNeurons(imname, xl, yl, windowSize); % Automatic detection

%% Exclude objects outside of crop

% This section makes sure only objects in the cropped area are being
% counted for validation purposes.

R = vertcat(stats.Centroid);
Rx = (R(:,1)>xl).*(R(:,1)<xl+windowSize);
Ry = (R(:,2)>yl).*(R(:,2)<yl+windowSize);
Rz = logical(Rx.*Ry);
bBox = vertcat(stats(Rz).BoundingBox);

stats = stats(Rz);

clear R
clear Rx
clear Ry
clear Rz
clear bBox
%% Parameters for exclusions

% Size (area) exclusion
areas = [stats.Area];
sortedAreas = sort(areas);
nAreas = round(0.1*length(areas)); % Not in use.
minArea = sortedAreas(nAreas);
minArea = 1;
% This size exclusion is currently dysfunctional by setting the area to 1.
% This only excludes objects with pixel are below 1; therefore excludes
% nothing. A functioning size filter is located in a dependency.

BoxRatio = 100.5;
% A parameter that excludes objects based on their bounding box ratio in an
% attempt to separate clusters from single cell objects. With current set
% of images this is not useful, so the ratio is set very high in order to
% make it non-functional.

singleCellArea = 130;
electrodeArea = 330;
% Approximate pixel areas of a single cell soma and a single electrode.

%% Plot detected objects on original image

% Plot the original image
imshow(I);
hold on;

% Plot detected regions
ctr = 0;
for idx = 1:length(areas)
    
    % Plot single neurons
    if stats(idx).Area <= singleCellArea && stats(idx).Area > minArea % only plot regions that satisfy size constraints (currently excludes nothing)
        box = stats(idx).BoundingBox(3:4);
        if ~(box(1) > box(2) * BoxRatio || box(2) > box(1) * BoxRatio) % only plot regions that satisfy shape constraints (all objects currently satisfy this)
            h = rectangle('Position', stats(idx).BoundingBox, 'EdgeColor', [0.75 0 0]); % Creates box around object.
            hold on;
            ctr = ctr+1; % Adds one count to total count of neurons.
        end
        
    % Plot clusters
    elseif stats(idx).Area > singleCellArea % Objects that are bigger than a single cell area; considered clusters.
        rectangleArea = stats(idx).BoundingBox(3) * stats(idx).BoundingBox(4);
        h = rectangle('Position',stats(idx).BoundingBox);
        set(h,'EdgeColor',[0.5 0 0], 'linewidth', 1.4);
        hold on;
        ctr = ctr + (stats(idx).Area/singleCellArea); % Adds average number of cells in a cluster by dividing the cluster pixel area by the average cell size.)
     end
end
% ctr = ctr - 59*electrodeArea/singleCellArea; % This line normally
% excludes electrode area from the count, however since crop selection is
% random, we don't know how many electrodes will be in an image. Electrode
% area has to be deducted manually during validation.
title(['There are ',  num2str(round(ctr)),  ' objects in this image']);
hold on;

% Plot the manually marked cells in green squares. These are all objects we
% clicked on initially that we think are neuronal soma.
scatter(ROI(:,1),ROI(:,2), 100, 's', 'markeredgecolor', [0 0.5 0], 'linewidth', 1.5);
hold on

% Plot regions that were correctly detected in yellow triangles. This is
% agreement between our manual marking and the scripts marking.
correct = parseCoord(ROI, stats, xl, yl, windowSize);
scatter(correct(:,1),correct(:,2), 20, 'y^', 'filled');

xlim([xl xl+windowSize])
ylim([yl yl+windowSize])
% Resets view of image to the original crop we are currently analysing.

