function correct = parseCoord(ROI, stats, xl, yl, windowSize)

% Parse the manual & automatic object detection results to get an estimate
% of the detection algorithm. 

% INPUT:
%   ROI - matrix containing the coordinates of the manually marked regions
%   stats - structure containing the object detection results
%   xl - scalar specifying the X coordinate of the cropped image
%   yl -  scalar specifying the Y coordinate of the cropped image
%   windowSize - scalar specifying the size of the image to be analyzed

% OUTPUT:
%   correct - coordinates of the correctly detected regions of interest

ROImanual = sort(ROI,1);

% Get only objects within the cropped image
R = vertcat(stats.Centroid);
Rx = (R(:,1)>xl).*(R(:,1)<xl+windowSize);
Ry = (R(:,2)>yl).*(R(:,2)<yl+windowSize);
Rz = logical(Rx.*Ry);
bBox = vertcat(stats(Rz).BoundingBox);

correct = [];

for i = 1:length(ROI)
    roi_cords = ROI(i,:);
    for j = 1:length(bBox)
        rectangle_cords = bBox(j,:);
        if roi_cords(1) > rectangle_cords(1) && roi_cords(1) < rectangle_cords(1)+rectangle_cords(3) && roi_cords(2) < rectangle_cords(2)+rectangle_cords(4) && roi_cords(2) > rectangle_cords(2)
            correct(end+1,:) = roi_cords;
        end
    end
end

[~, icorrect] = uniquetol(correct, 'byrows', true);
correct = correct(sort(icorrect), :);

end