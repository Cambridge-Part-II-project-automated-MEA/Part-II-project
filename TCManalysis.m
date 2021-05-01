
%% Total Cell Matter Area calculator

% Calculates all area of an image that is covered in cell matter. Currently
% not in use, but could be used anytime.

files = dir('*.jpg');

tic
for file = 1:length(files)
    I = imread(files(file).name);
    % Igrey = rgb2gray(I);
    Icomp = imcomplement(I);
    BW = imbinarize(Icomp,'adaptive', 'ForegroundPolarity','bright','Sensitivity', 0.55);
    stats1 = regionprops(BW, 'Area');
    TCMarea(file) = sum([stats1.Area]);
    PercentageTCM(file) = TCMarea(file)/(width(BW)*length(BW));
end
toc

filenames = {files.name};

results = table(filenames', TCMarea', PercentageTCM', 'VariableNames', {'culture';'TCMarea'; 'PercentageTCM'});

save('results.mat', 'filenames', 'results');

strfind(results{:,1},'richn10')

pos = strfind(results{:, 1},'richn10');

cell2mat(pos)

% results{1, :}

idx = find(contains(results.culture,'richn7'));
