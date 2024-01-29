function [cellMap,nROI] = cellDetectFun(tiffData,minSiz,Tbin)
% This function takes a TIF stack in input and detects the active ROIs
% (such as astrocytes or neurons) present over the entire recording.
% ----------------------------------------------------------------------- %
% *** Inputs ***
% * tiffData * = TIF file containing the fluorescent imaging data
% * minSiz *   = minimal size (in pixels) of a detected cell
% * Tbin *     = adaptive binarization threshold parameter (default 0.5) 
% ----------------------------------------------------------------------- %
% *** Outputs ***
% * cellMap * is an image of the same size as the frames of *tiffData*
% giving the position of the detected ROIs in the image and associating 
% them with their label number
% * nROI * is the total number of detected ROIs
% ----------------------------------------------------------------------- %
% L. Zonca, Jan. 2022
% ----------------------------------------------------------------------- %

% Sum all frames over time
sumData = sum(tiffData,3);
% Normalize the summed data
normData = (sumData-min(sumData(:)))./(max(sumData(:))-min(sumData(:)));
% Detect all astros
if ~exist('Tbin','var') || isempty(Tbin)
    Tbin = 0.5;
end
T_full = adaptthresh(normData,Tbin);
astros_full = imbinarize(normData,T_full);

% Label all astros
[cellMap,nROI] = label_2photons(astros_full,minSiz);

end