function GloEvMat = globalEvtActiv(GlobalEventsTime,realEvts,varList)
% This function builds a cell array containing the list of active ROIs and
% their peak time for each global event based on *GlobalEventsTime* from
% the function globalEvSeg
% ----------------------------------------------------------------------- %
% *** Inputs ***
% * GlobalEventsTime * is a vector containing the time points (in frames) of 
% the average signal at the detected local minima 
% * realEvts * is the cell array conatining all event features for each active 
% ROI. It is the output of the function segWithBaseCorr
% * varList * is the list of active ROIs, it is an output of the function segWithBaseCorr
% ----------------------------------------------------------------------- %
% *** Outputs ***
% * GloEvMat * is a cell array containing the list of active ROIs and
% their peak time for each global event. It is used to build the activation
% paths and the connectivity graphs.
% ----------------------------------------------------------------------- %
% L. Zonca, Jan. 2022
% ----------------------------------------------------------------------- %

GloEvMat = cell(length(GlobalEventsTime),2);
for ev = 1:length(GlobalEventsTime)
    for a = 1:length(varList)
        [isActive,locEv] = ismember(ev,realEvts{a,4});
        if isActive == 1
            GloEvMat{ev,1} = [GloEvMat{ev,1}; varList(a)]; % List of active astrocytes in event
            GloEvMat{ev,2} = [GloEvMat{ev,2}; realEvts{a,2}(locEv)]; % Time of peak for each astrocyte in event
        end
    end
end
end
