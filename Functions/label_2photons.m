function [labeled,nROI] = label_2photons(detected,minSiz)
% This function labels the detect ROIs (cells) and remove the ROIs smaller
% than minSiz pixels it is called in the cellDetectFun function to label.
% ----------------------------------------------------------------------- %
% *** Inputs ***
% * detected * is a binary image giving the ROIs to filter and label
% * minSiz * is the minimal size (in pixels) of a ROI that will be kept 
% ----------------------------------------------------------------------- %
% *** Outputs ***
% * labeled * is an image of the same size as *detected* giving the position 
% of the ROIs in the image and associating them with their label number
% * nROI * is the total number of detected ROIs
% ----------------------------------------------------------------------- %
% L. Zonca, Jan. 2022
% ----------------------------------------------------------------------- %

%%%%%%%%% DETECTION OF CONNECTED COMPONENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = bwconncomp(detected,8);       
labeled = labelmatrix(CC);

%%%%%%%%% FILTERING TOO SMALL ZONES (<sizMin sq. pix) %%%%%%%%%%%%%%%%%%%%%%%%%
tailleAstro = cellfun(@numel, CC.PixelIdxList);
for k=1:length(CC.PixelIdxList);
    if tailleAstro(k)<minSiz
       labeled(CC.PixelIdxList{k})=0;
       CC.PixelIdxList{k}=[];
    end
end

idx = cellfun('isempty',CC.PixelIdxList);
CC.PixelIdxList(idx) =  [];
CC.NumObjects = length(CC.PixelIdxList);

labeled = labelmatrix(CC);
%%%%%%%%% DISPLAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
astro_label = label2rgb(labeled,'jet', 'w', 'shuffle');
figure; 
imshow(astro_label)

%%% Display of the label numbers on the centroid of each object %%%
s = regionprops(labeled, {'Centroid', 'PixelIdxList'});
% Make the image uint8
I = im2uint8(astro_label);
I(~astro_label)=200;
imshow(I, 'InitialMagnification', 'fit')
% Plot the number of each object at the corresponding centroid
hold on
for k = 1:numel(s)
   centroid = s(k).Centroid;
   text(centroid(1), centroid(2), sprintf('%d', k));
end
hold off
% Number of detected cells
nROI = CC.NumObjects; 
end