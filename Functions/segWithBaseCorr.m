function [meanFreq,realEvts,diffToBaseAll,varList] = ...
    segWithBaseCorr(signalMat,dt,dtev,minDepth,maxEvThresh,Tev,cfDown,smoothParam,evNumVec,GlobalEventsTime)
% This function segments the input signals in signalMat and extract the
% segmented events and inter-events features detailed in *Outputs*.
% ----------------------------------------------------------------------- %
% *** Inputs ***
% * signalMat *   = matrix containning the signals to segment
% each trace must be in one column i.e. size(signalMat) = N x nT where N is
% the number of time samples and nT the number of traces.
% with N the length in time and nT the number of traces
% * dt *          = signal time sampling
% * dtev *        = caracteristic event time 
% (to be removed around peaks to build the signal outside of peaks)
% * minDepth *    = min variation level to detect the peaks
% * maxEvThresh * = min threshold for event detection if std(noise) is too high
% * Tev *   = coef to noise for event detection threshold
% then the event detection threshold is t_ev = min(noiseCoef*std(signal noise),minEvThresh)
% * cfDown *      = coef to compensate slower decay than increase for the events
% (to build the signal outside of peaks)
% * smoothParam * = smoothing level for the spline for baseline computation
% * evNumVec / GlobalEventsTime * OPTIONAL = Global events number (outputs
% from globalEvSeg) - NECESSARY for Graphs and spatial correlations
% analysis
% ----------------------------------------------------------------------- %
% *** Outputs ***
% meanFreq      = mean event frequency for each trace in signalMat
% realEvts{:,1} = beginning of events
% realEvts{:,2} = peak of events
% realEvts{:,3} = end of events
% realEvts{:,4} = global event number (for spatial analysis in image data)
% realEvts{:,5} = events amplitude
% realEvts{:,6} = number of sub peaks
% realEvts{:,7} = sub peak frequency
% varList       = list of traces with active events (vector)
% diffToBaseAll = matrix of size(*signalMat*) containing the corrected
% signals (i.e. with baseline correction)
% ----------------------------------------------------------------------- %
% L. Zonca, Jan. 2022
% ----------------------------------------------------------------------- %

[N,nT] = size(signalMat);
% Colormap for the plots
cmap = hsv(max(10,nT)); 

% Time vector
timeVec = 1:N;
% Global envents info (for spatial analysis when image data)
if ~exist('evNumVec','var') || ~exist('GlobalEventsTime','var') 
    evNumVec = zeros(size(timeVec));
    GlobalEventsTime = [1 N];
else
    GlobalEventsTime = GlobalEventsTime';
end

% Initiate outputs
meanFreq = zeros(nT,1);
realEvts = cell(nT,7);
varList = 1:nT;
diffToBaseAll = zeros(size(signalMat));

for traceIdx = 1:nT
    signal = signalMat(:,traceIdx);
    % plot the trace
    figure
    plot(timeVec*dt,signal,'LineWidth',1,'color',cmap(traceIdx,:))
    % Detect peaks
    noEvtsPeriods = timeVec;
    [eventsVal,eventsTime] = findpeaks(signal,'MinPeakProminence',minDepth);
    % Extract signal outside of peaks for baseline fit
    for ev = 1:length(eventsTime)
        if eventsTime(ev)>dtev && eventsTime(ev)+cfDown*dtev<=timeVec(end)
            noEvtsPeriods = setdiff(noEvtsPeriods, eventsTime(ev)-dtev:eventsTime(ev)+cfDown*dtev);
        elseif eventsTime(ev)>dtev
            noEvtsPeriods = setdiff(noEvtsPeriods, eventsTime(ev)-dtev:timeVec(end));
        elseif eventsTime(ev)+cfDown*dtev<=timeVec(end)
            noEvtsPeriods = setdiff(noEvtsPeriods, 1:eventsTime(ev)+cfDown*dtev);
        end
    end
    % Process only if detected events
    if ~isempty(eventsTime) 
        % Correct boundaries of the signal outside of peaks at beginning and end
        if eventsTime(1)<dtev
            noEvtsPeriods = [1 noEvtsPeriods];
        end
        noEvtSig = zeros(size(signal));
        noEvtSig(noEvtsPeriods,:) = signal(noEvtsPeriods,:);
        noEvtSig(GlobalEventsTime,:) = ...
            repmat((mean(signal)-min(signal))/2+min(signal),...
            length(GlobalEventsTime),1);
        if ~ismember(length(signal(:,1)),noEvtsPeriods)
            noEvtSig(end,:) = (mean(signal)-min(signal))/2+min(signal); 
            noEvtsPeriods = [union(noEvtsPeriods,GlobalEventsTime), length(signal(:,1))];
        else noEvtsPeriods = union(noEvtsPeriods,GlobalEventsTime);
        end
        % --- OUPUT 1 --- Mean oscillation frequency (between peaks) ---
        if length(eventsTime)>=2
            periodes = (eventsTime(2:end)-eventsTime(1:end-1))*dt;
            frequences = 1./periodes;
            meanFreq(traceIdx) = mean(frequences);
        end
        hold on
        % Plot detected peaks
        plot(eventsTime*dt,eventsVal,'.','MarkerSize',20,'color',cmap(traceIdx,:))
        % Plot signal used for baseline fit
        plot(noEvtsPeriods*dt,noEvtSig(noEvtsPeriods,:),'LineWidth',1.5,'color','r')
        % Fit the baseline outside of events
        baseLine = fit(noEvtsPeriods',noEvtSig(noEvtsPeriods,:),'smoothingspline','SmoothingParam',smoothParam);
        % Plot the fit
        plot(noEvtsPeriods*dt,baseLine(noEvtsPeriods),'k','LineWidth',1.5)
        baseLineInt = interp1(unique(noEvtsPeriods),baseLine(unique(noEvtsPeriods)),timeVec,'nearest');
        % Signal correction: remove baseline
        diffToBase = signal-baseLineInt';
        diffToBaseAll(:,traceIdx) = diffToBase;
        % Plot corrected signal
        plot(timeVec*dt,signal-baseLineInt','color','g','LineWidth',1.5)
        plot(timeVec*dt,zeros(size(timeVec)),'k','LineWidth',1.5)
        
        % --- OUTPUT 2 --- beginning-peak-end times / amplitude / number of supeaks - sub peak frequency ---
        % Detect event boundaries using corrected signal (diffToBase)
        signDiffToBase = diffToBase(2:end).*diffToBase(1:end-1);
        newEvents = find(signDiffToBase <= 0);

        realEvts{traceIdx,1} = []; % beginning of events
        realEvts{traceIdx,2} = []; % peak of events
        realEvts{traceIdx,3} = []; % end of events
        realEvts{traceIdx,4} = []; % global event number (for spatial analysis)
        realEvts{traceIdx,5} = []; % events amplitude
        realEvts{traceIdx,6} = []; % number of sub peaks inside events
        realEvts{traceIdx,7} = []; % sub peak frequency
        % Define minimal event amplitude proportional to noise std dev of this astrocyte
        noiseAmpli = std(diffToBase(noEvtsPeriods));
        minEvAmpli = min(Tev*noiseAmpli,maxEvThresh);
        for ev = 1:length(newEvents)
            [evVal,evTime] = max(diffToBase(newEvents(max(1,ev-1)):newEvents(ev)));
            if evVal > minEvAmpli % If peak amplitude is too small, discard event
                realEvts{traceIdx,1} = [realEvts{traceIdx,1}; newEvents(max(1,ev-1))];
                realEvts{traceIdx,2} = [realEvts{traceIdx,2}; newEvents(max(1,ev-1))+evTime-1];
                realEvts{traceIdx,3} = [realEvts{traceIdx,3}; newEvents(ev)];
                realEvts{traceIdx,4} = [realEvts{traceIdx,4}; evNumVec(newEvents(max(1,ev-1))+evTime-1)];
                realEvts{traceIdx,5} = [realEvts{traceIdx,5}; diffToBase(newEvents(max(1,ev-1))+evTime-1)];
                try 
                    [subPeakVal,subPeakTime] = findpeaks(diffToBase(newEvents(max(1,ev-1)):newEvents(ev)),...
                        'MinPeakProminence',minEvAmpli/1.5);
                    if ~isempty(subPeakTime)
                        realEvts{traceIdx,6} = [realEvts{traceIdx,6}; numel(subPeakVal)];
                        plot(timeVec(newEvents(max(1,ev-1))+subPeakTime-1)*dt,...
                            diffToBase(newEvents(max(1,ev-1))+subPeakTime-1),...
                        'x','color',cmap(mod(ev,length(cmap)),:),'MarkerSize',8,'LineWidth',2)
                        if numel(subPeakVal)>1
                            realEvts{traceIdx,7} = [realEvts{traceIdx,7}; 1/(mean((subPeakTime(2:end)-subPeakTime(1:end-1)))*dt)];
                        end
                    end
                catch
                    warning(['Skipped evt ' num2str(ev)])
                end
            end
        end 
        plot(realEvts{traceIdx,2}*dt,realEvts{traceIdx,5}+baseLineInt(realEvts{traceIdx,2})','.','MarkerSize',20,'color',cmap(traceIdx,:))
        title(num2str(minEvAmpli))
        plot(noEvtsPeriods*dt,baseLine(noEvtsPeriods)+minEvAmpli,'k--','LineWidth',1.)
        plot(noEvtsPeriods*dt,baseLine(noEvtsPeriods)-minEvAmpli,'k--','LineWidth',1.)
        xlim([min(timeVec*dt),max(timeVec*dt)])
        xlabel('Time (s)')
        ylabel(['Trace # ' num2str(traceIdx)])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Keep this trace only if it has active events
        if length(realEvts{traceIdx,2})<1
            varList = setdiff(varList,traceIdx);
            close
        end
    else varList = setdiff(varList,traceIdx);
         close
    end
end