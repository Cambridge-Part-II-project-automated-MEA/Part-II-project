
% Count total area of binarized image: Total Cell Matter

stats1 = regionprops(BW, 'Area'); % Lists all connected components including one pixel
% sized objects, meaning all black pixels are included.
% This uses an already binarized image, which uses a chosen parameter of
% sensitivity. This needs to be adjusted according to needs or an objective
% metric, such as weighted darkness of all pixels.
TCMarea = sum([stats1.Area]);
% Sums the size of all above objects, giving total count of pixels.
PercentageTCM = TCMarea/(width(BW)*length(BW));
% Total percentage area of photograph taken up by Cell Matter: achieved by
% dividing total pixels by the product of the resolution of the image.

%% Superimposing binarized image onto original

Iorg = imread('richnhighres1.jpg');
% Image to be analysed

% Igrey = rgb2gray(I); % Use if input image is not already greyscale, otherwise comment out

% Icomp = imcomplement(Igrey); % Change Igrey to I if not using rgb2grey

BW = imbinarize(Icomp,'adaptive', 'ForegroundPolarity','bright','Sensitivity', 0.55);
% The number between 0-1 is the sensitivity, higher number means more
% pixels counted, the area it considers cell matter area increases

% se = strel('disk', 2); % Use if small objects need to be excluded.
% Smaller number --> see more
% The higher the number, the bigger objects it will exclude

% Iopenned = imopen(BW,se); % Use if strel functio is used

figure, imshowpair(BW, Iorg); % Change BW to Iopenned if strel is used.

%% 

% Count total area of binarized image: Total Cell Matter

stats1 = regionprops(BW, 'Area'); % Lists all connected components including one pixel
% sized objects, meaning all black pixels are included.
% This uses an already binarized image, which uses a chosen parameter of
% sensitivity. This needs to be adjusted according to needs or an objective
% metric, such as weighted darkness of all pixels.
TCMarea = sum([stats1.Area]);
% Sums the size of all above objects, giving total count of pixels.
PercentageTCM = TCMarea/(width(BW)*length(BW));
% Total percentage area of photograph taken up by Cell Matter: achieved by
% dividing total pixels by the product of the resolution of the image.

