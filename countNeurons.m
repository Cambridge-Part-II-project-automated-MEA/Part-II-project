function stats = countNeurons(imname, xl, yl, windowSize)
%% Counting Neurons'
% 
% INPUT:
%   imname - string specifying the image to be analysed
%   xl - scalar specifying the X coordinate of the cropped image
%   yl -  scalar specifying the Y coordinate of the cropped image
%   windowSize - scalar specifying the size of the image to be analyzed

% OUTPUT:
%   stats - structure containing the object detection results

%% Read Image
I = imread(imname);

% If not cropping the image, set default values to full image
if ~exist('xl','var')
    xl = size(I, 1);
end
if ~exist('yl','var')
    yl = size(I, 2);
end
imshow(I);
title('Select 6 consecutive electrodes');

for i = 1:6 % Click to select six electrodes.
    ROI = drawpoint;
    TemplateElec(i,:) = ROI.Position;
    clear ROI
end

for i = 1:length(TemplateElec) - 1
    lengthx(i) = TemplateElec(i+1,1)-TemplateElec(i,1);
end

eleclength = mean(lengthx); % Calculates mean of the length between electrodes.

Elec1= TemplateElec(1,:);
Elec1(1)= Elec1(1)- eleclength;

count = 0;
for w = 1:8  %w=width cooresponding to y coord
    for l = 1:8 %l=length corresponding to x coord
        if (l==1 && w==1) || (l==8 && w==8) || (l==1 && w==8) || (l==8 && w==1) || (l==5 && w==1)
            continue
        else
            count=count+1;
            ElecCoord(count,1)= Elec1(1)+ (w-1)*eleclength;
            ElecCoord(count,2)= Elec1(2)+ (l-1)*eleclength;
        end
    end
end

ElecCoord = ElecCoord';

% This whole segment calculates location of electrodes based on an initial
% 6 clicks. This is not necessary for calculation of number of cells on an
% image.
%% Image enhancement
I = imread(imname);
% Reads in image.
ref = imread('richn1.jpg'); % Reference image with high contrast, which will be used 
% for histogram matching.

I = imsharpen(I, 'threshold', 0.2, 'amount', 2, 'radius', 4); 
% Sharpens image: parameters can be changed to match need.
I = imhistmatch(I, ref, 'method', 'uniform');
% Histogram matching: matches the pixel shade distribution of an image
% to that of a reference image provided previously.
Icomp = imcomplement(I);
% Inverts the colour of the image.
Icomp = imsharpen(Icomp, 'threshold', 0.1, 'amount', 2, 'radius', 1);
% Another round of sharpening the image.

%%

WeightedDarknessScript; 
% Calculates weighted darkness index of greyscale picture.
% Currently not in use. Might later be used as an automatic way of
% adjusting sensitivity of binarization.

ScaledSensitivityBinarization = WeightedDarkness*0.86;
% Currently not in use: an attempt to calculate a sensitivity based on
% overall shade of the image.

BW = imbinarize(Icomp,'adaptive', 'ForegroundPolarity','bright','Sensitivity', 0.05);
% Binarizes the image. Sensitivity can be adjusted to match need, but
% current script has been validated based on the 0.05 value.


se = strel('disk', 2); 
% Creates a mask that will filter out small objects such as debris.
% Change number after 'disk' to a smaller: leaves smaller objects in image.
% Larger number: makes more objects of increasing size disappear.
Iopenned = imopen(BW,se);
% Applies previously made mask.
figure, imshowpair(Iopenned, I);
imshow(Iopenned);
% Shows current image progress.

CC = bwconncomp(Iopenned, 4);
% Counts connected components; 'creates' objects from connected 'islands'
% of pixels.
stats = regionprops(CC, 'Eccentricity', 'Area', 'BoundingBox', 'Image', 'Centroid');
% Extracts information about the objects in the image and collects them in
% a table.
stats = stats([stats.Area] > 21); 
% Removes all objects below given pixel Area size. This filters out much
% debris. Change number to match need if using a different type of culture.
areas = [stats.Area];
% Creates list of remaining objects that weren't removed due to small size.
numObjects = length([stats.Area]);
% Counts above objects.
BoundingBox = [stats.BoundingBox];
for i = 1:length(areas)
    BoundBox(i,:) = BoundingBox (1 + (i * 4 -4): i*4);
end
clear BoundingBox
TotalCellMatterArea;
% Calculates percentage of image area covered in cell matter.
AreaValsBig = [stats.Area];
colsToDelete = AreaValsBig < 110;
AreaValsBig(colsToDelete) = [];
% Creates list of objects that are bigger than given value (Approximate
% pixel size of a single cell soma).
AvgCellsPerCluster = AreaValsBig/110;
% Approximates number of cells in each cluster based on their pixel size.
TotalClusterCells = sum(AvgCellsPerCluster);
% Sums number of cells in all clusters.
ElecArea = 59*330/110;
% Subtracts a constant number from above count that represent the
% electrodes, since the 59 electrodes are all counted as clusters too.
% 59 electrodes, each 330 pixels in size, divided by the average pixel size
% of a single cell soma.
FinalnumObjects = numObjects - length(AreaValsBig) + TotalClusterCells - ElecArea;
% Calculates final number of objects.
end
