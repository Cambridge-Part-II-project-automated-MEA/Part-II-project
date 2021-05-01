%% Batch analysis
% This script is simplified for quicker analysis. None of the analysed
% images are shown visually, this batch analysis script only calculates and
% collects the information about number of cells in the image.


% Dependencies:
% Reference image

clearvars
% Clears all variables.

files = dir('*.jpg');

% Centroids = cell(5000, 5000);
% Centroids(:,:) = {0};

tic % Will measure time taken to run script.
for file = 1:length(files) % Sets length of for loop based on number of images to be analysed.
    imname = (files(file).name);
    
    I = imread(imname);
    ref = imread('richn1.jpg'); % Reference image for histogram matching.
    I = imsharpen(I, 'threshold', 0.2, 'amount', 2, 'radius', 4); % Sharpening of image.
    I = imhistmatch(I, ref, 'method', 'uniform'); % Matching the pixel shade distribution histograms to reference image.
    Icomp = imcomplement(I); % Inverts image shade.
    Icomp = imsharpen(Icomp, 'threshold', 0.1, 'amount', 2, 'radius', 1); % Sharpens image.
    BW = imbinarize(Icomp,'adaptive', 'ForegroundPolarity','bright','Sensitivity', 0.05); % Binarizes image making it black and white.
    se = strel('disk', 2); % Creates mask that removes small objects such as debris.
    Iopenned = imopen(BW,se); % Applies the mask.
    CC = bwconncomp(Iopenned, 4); % Counts connected components; 'creates' the objects so to say.
    stats = regionprops(CC, 'Area', 'Centroid'); % Creates table with information about the objects detected.
    stats = stats([stats.Area] > 21); % Removes all objects below given pixel Area size.
    areas = [stats.Area]; % Creates new list of objects without the small objects removed in previous step.
    numObjects = length([stats.Area]); % Counts number of objects.
    AreaValsBig = [stats.Area];
    colsToDelete = AreaValsBig < 110;
    AreaValsBig(colsToDelete) = []; % Creates list with only clusters in it (objects bigger than approximate size of a single cell soma)
    AvgCellsPerCluster = AreaValsBig/110; % Approximates nunmber of cells in each cluster.
    TotalClusterCells = sum(AvgCellsPerCluster); % Sums all cell counts of clusters.
    ElecArea = 59*330/110; % Constant representing the approximate number of cells the electrodes add.
    FinalnumObjects(file) = round (numObjects - length(AreaValsBig) + TotalClusterCells - ElecArea); % Final cell count after removing electrodes.
%     Centroids(:, file) = [stats.Centroid];
%     Centroids(:, file) = {stats.Centroid}.';
    % Currently the script does not need to save the location of each
    % neuron on each image, therefore commented out and not yet complete.
end
toc % Measures time needed for script to batch analyse all images.

filenames = {files.name};

results = table(filenames', FinalnumObjects', 'VariableNames', {'culture';'Number of cells'});

save('results.mat', 'filenames', 'results');
% Creates table of results.

% strfind(results{:,1},'richn10')

% pos = strfind(results{:, 1},'richn10');

% cell2mat(pos)

% results{1, :}

idx = find(contains(results.culture,'richn7'));
% Used to find specific image information.

