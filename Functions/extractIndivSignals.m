function [decharges,dechNorm,varList] = extractIndivSignals(tiffData,cellMap,dSiz,ROI_list,varMin)
% This function extracts the fluorescence signal of each ROI in *ROI_list*
% over the entire recording *tiffData* and removes the signals with less than
% *varMin* amplitude variations in their normalized signal compared to
% their min value
% ----------------------------------------------------------------------- %
% *** Inputs ***
% * tiffData * = TIF file containing the fluorescent imaging data
% * cellMap *  = image of the size of one frame of *tiffData* containing
% the labeled ROIs from which you wish to extract the signals. The ROI
% labels should be consistent with those in *ROI_list*
% * dSiz *     = spatial down sampling factor for the extraction (dSiz = 1
% means no down sampling, dSiz should not be larger than the minimal ROI
% size)
% * ROI_list * = vector containing the list of the ROI labels from which
% you wish to extract the signal
% * varMin *   = number between 0 and 1 representing the minimal fluctation
% level of the signal (compared to its mean value) to consider the ROI
% as active
% ----------------------------------------------------------------------- %
% *** Outputs ***
% * decharges * = Matrix of size length(*tiffData*) (i.e. total duration of
% recording) x length(*ROI_list*) (i.e. number of considered ROIs) containing
% the signal extracted for each ROI
% * dechNorm *  = Matrix of size length(*tiffData*) (i.e. total duration of
% recording) x length(*ROI_list*) (i.e. number of considered ROIs) containing
% the nomralized signal (between 0 and 1) extracted for each ROI
% * varList *  = vector containing the labels of the active ROIs (i.e. with
% variations greater than *meanVar*)
% ----------------------------------------------------------------------- %
% L. Zonca, Jan. 2022
% ----------------------------------------------------------------------- %

% signal extraction
D = 1;                        % beginning
F = length(tiffData(1,1,:));  % end of signal
decharges = zeros(F-D+1,length(ROI_list));

for A = 1:length(ROI_list)
    [x,y] = find(cellMap == ROI_list(A));
    x_ssEch=min(x):dSiz:max(x);
    y_ssEch=min(y):dSiz:max(y);
    surface=length(x_ssEch)*length(y_ssEch);
    for j = D:F-D+1
       % Take mean luminosity over each astrocyte's surface
       decharges(j,A)=sum(sum(tiffData(x_ssEch,y_ssEch,j)))/surface; 
    end
end

% Normalization of the extracted signals and filtering of inactive ROIs
dechNorm = zeros(size(decharges));
varList = [];
for nas = ROI_list
    % Normalize signal
    dechNorm(:,nas) = (decharges(:,nas)-min(decharges(:)))/(max(decharges(:))-min(decharges(:)));
    % Plot only if variations/mean value >= varMin
    var = max(dechNorm(:,nas)-mean(dechNorm(:,nas)));
    if var >= varMin
        % Keep astro number in list of active astrocytes
        varList = [varList; nas];
    end
end
end