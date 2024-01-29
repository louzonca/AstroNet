function [cr,corrCurve,corThresh] = pairCrossCorr(diffToBaseAll,varList,dT)
% This function builds the correlation curve based on the pairwise Pearson
% correlation coefficient (computed by the Matlab function corr) of all the
% corrected signals in *diffToBaseAll*
% ----------------------------------------------------------------------- %
% *** Inputs ***
% * diffToBaseAll * is a matrix of (size length of time recordings) x (number 
% of traces) containing signals with baseline correction (from function
% segWithBaseCorr)
% * varList * is the list of active ROIs
% * dT * is a value between 0 and 1 giving the step size used to build the
% correlation curve
% ----------------------------------------------------------------------- %
% *** Outputs ***
% * cr * is the correlation matrix of diffToBaseAll
% * corrCurve * is a vector containing the proportion of signals with at 
% least one correlation over each corresponding value in *corThresh* i.e 
% plot(corThresh,corCurve) plots the correlation curve
% ----------------------------------------------------------------------- %
% L. Zonca, Jan. 2022
% ----------------------------------------------------------------------- %

cr = corr(diffToBaseAll);
% Obtain correlation curve
corThresh = 0:dT:1;
corrCurve = zeros(size(corThresh));
for cTidx = 1:length(corThresh)
    cT = corThresh(cTidx);
    [a,~]=ind2sub(size(cr),find(cr-eye(size(cr))>cT));
    corrCurve(cTidx) = numel(unique(a))./length(varList);
end
end