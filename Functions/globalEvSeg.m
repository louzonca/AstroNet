function [GlobalEventsVal,GlobalEventsTime,evNumVec] = globalEvSeg(dechNorm,varList,filtWndw,minDepth)
% This function detects the local minima in the average signal (over all
% ROIs) to segment it into 'global' events
% ----------------------------------------------------------------------- %
% *** Inputs ***
% * dechNorm * is a matrix containing the normalized (between 0 and 1)
% signals of all considered ROIs. *dechNorm* can be obtained from 
% the extractIndivSignals function
% * varList * is the vector containing the labels of the active ROIs it can
% be obtained from the extractIndivSignals function
% * filtWndw * is the size (in frames) of the window used to low-pass
% filter the averaged signal (over all ROIs) (larger values lead to more
% filtering)
% * minDepth * is the minimal depth for a local minima (used to segment
% into 'global' events)
% ----------------------------------------------------------------------- %
% *** Outputs ***
% * GlobalEventsVal * is a vector containing the value of the average signal 
% at the detected local minima
% * GlobalEventsTime * is a vector containing the time points (in frames) of 
% the average signal at the detected local minima 
% * evNumVec * is a vector of the length of the averaged signal containing
% the number of the global event at each time
% ----------------------------------------------------------------------- %
% L. Zonca, Jan. 2022
% ----------------------------------------------------------------------- %

% Mean over all active astrocytes
meanSig = mean(dechNorm(:,varList),2);
meanSigFilt = meanSig;
% Low pass filter
for i = filtWndw/2+1:length(meanSig)-filtWndw/2
    meanSigFilt(i) = mean(meanSig(i-filtWndw/2:i+filtWndw/2));
end

% Find local minima to segment the global events
[GlobalEventsVal,GlobalEventsTime] = findpeaks(-meanSigFilt,'MinPeakProminence',minDepth);

timeVec = 1:length(dechNorm(:,1));
evNumVec = zeros(size(timeVec));
for e = 1:length(GlobalEventsTime)-1
    evNumVec(GlobalEventsTime(e):GlobalEventsTime(e+1))= e;
end
evNumVec(GlobalEventsTime(end):end)=length(GlobalEventsTime);
end