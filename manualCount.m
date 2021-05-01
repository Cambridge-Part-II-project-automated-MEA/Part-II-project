function [ROI, xl, yl] = manualCount(imname, windowSize)

% This script serves to facilitate the manual counting of neurons
%  If we save the variables: xl, yl, ROI
%  We will be able to improve the benchmarking of the detection script to
%  include not only the number of objects/area but also the exact
%  sensitivity and precision. This script was used to validate the
%  detection script.

% INPUT:
%   imname - string specifying the image to be analysed.
%   windowSize - scalar specifying the size of the image to be analyzed.

% OUTPUT:
%   ROI - matrix containing the coordinates of the manually marked regions.
%   xl - scalar specifying the X coordinate of the cropped image.
%   yl -  scalar specifying the Y coordinate of the cropped image.

% Note: [xl, xl+windowSize; yl, yl+windowSize] are the full coordinates of
% the cropped image

I = imread(imname);
% Reads in image.
imshow(I);
title('Select regions of interest and click "Done"')

xl = randi([1,size(I, 1)-windowSize],1);
yl = randi([1,size(I, 2)-windowSize],1);
% Creates random coordinates to create a cropped part of the image.

xlim([xl xl+windowSize])
ylim([yl yl+windowSize])

button = uicontrol('Position',[300 15 150 30],'String','Done',...
    'Callback', @Done);
% Closes image once finished selecting neuronal soma.

flg = 0;
counter = 1;
while ~flg
    roi = drawpoint('color',[0 0.5 0]);
    if ~isvalid(roi) || isempty(roi.Position)
        flg = 1;
    else
        ROI(counter, :) = roi.Position;
        counter = counter+1;
    end
end
    function Done(button, EventData)
        close(gcf);
    end
end
% Counter. Clicking with the mouse on the image on areas which we think are
% neuronal soma. This will save those co-ordinates for later comparison
% with how the script performs. With this comparison we can calculate the
% precision and sensitivity of the script.
