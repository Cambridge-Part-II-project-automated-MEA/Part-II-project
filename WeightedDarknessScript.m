
% Weighted Darkness of Greyscale Picture

SumColumns = sum(Icomp);
SumImage = sum(SumColumns);
WeightedDarkness = SumImage/(255*(width(Icomp)*length(Icomp)));

% This script calculates the overall shade of the image.
% Output is a number between 0-1. (between all pixel values being 255 or
% all values being 0)
% Currently this metric is not in use.