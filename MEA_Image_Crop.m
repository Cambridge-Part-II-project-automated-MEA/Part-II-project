%% Crop image around electrodes
% Courtesy of I. Loureiro.


clearvars
close all
% Clear variables from workspace.
grayImage = imread('richnhighres10.jpg');
% Image name to be cropped.
% subplot(1,2,1);
p = [10 10 60 30];
    set(0, 'DefaultFigurePosition', p)
    f1 = figure;
imshow(grayImage);

for i = 1:6 % to select six electrodes
ROI = drawpoint;
TemplateElec(i,:) = ROI.Position; 
clear ROI
end

for i=1:length(TemplateElec) - 1
    lengthx(i) = TemplateElec(i+1,1)-TemplateElec(i,1);
end  %calculates the length of x for all electrodes chosen
 
eleclength=mean(lengthx); %calculates mean of the length between elecs

Elec1= TemplateElec(1,:);
Elec1(1)= Elec1(1)- eleclength;

count= 0;
for w = 1:8  %w=width cooresponding to y coord           
    for l = 1:8 %l=length corresponding to x coord
        if (l==1 && w==1) || (l==8 && w==8) || (l==1 && w==8) || (l==8 && w==1) 
            count=count+1;
            ElecCoord(count,1)= Elec1(1)+ (w-1)*eleclength;
            ElecCoord(count,2)= Elec1(2)+ (l-1)*eleclength; % gives corner coords of grid
        end
    end
end

ElecCoord = ElecCoord';

I2 = imcrop(grayImage,[ElecCoord(1,1)- eleclength/2, ElecCoord(2,1)-eleclength/2, ...
    (ElecCoord(1,3)-ElecCoord(1,1))+eleclength, ...
    (ElecCoord(2,2)-ElecCoord(2,1))+eleclength]);
% Crops image based on electrode locations.
f2 = figure;
imshow(I2)
% J = imresize(I2, 1.5);
% f3 = figure;
% imshow(I2)


% figname = 'bad';
% saveas(f3,strcat(figname,'.jpg'))
