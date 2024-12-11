% This script concatenates all the processed data contained in the
% allCTX and allCA1 mat files in Matlab arrays to make the plots of the
% individual statistics of the events (Fig. S1 in the manuscript) and the
% correlation curves of Fig. 3
% --- Load data --------------------------------------------------------- %
load('Data\allCTX.mat')
load('Data\allCA1.mat')

% Concatenate individual values for all sessions CTX ------------------------------- %
recList = 1:2:5;
% Initiate (CTX)
tUpCTX = [];
tDownCTX = [];
evtdurCTX = [];
evtFreqCTX = [];
evtAmpliCTX = [];
numSubEvtCTX = [];
freqSubEvtCTX = [];
interEvtCTX = [];
actNumCTX = [];
totNumCTX = [];

% Concatenate (CTX)
for recNum = recList
    dt = 0.5;
    for s = 1:length(allCTX(recNum,:))
        %%%% ALL CTRL SESSIONS CTX %%%%
        if ~isempty(allCTX{recNum,s})
            evtdurCTX = [evtdurCTX; (cell2mat(allCTX{recNum,s}.evtEnd) - cell2mat(allCTX{recNum,s}.evtBegin))*dt];
            for e = 1:length(allCTX{recNum,s}.evtBegin)
                if length(cell2mat(allCTX{recNum,s}.evtBegin(e)))>=2
                    beg = cell2mat(allCTX{recNum,s}.evtBegin(e));
                    nd = cell2mat(allCTX{recNum,s}.evtEnd(e));
                    interEvtCTX = [interEvtCTX; (beg(2:end)-nd(1:end-1))*dt];
                end
            end
            tUpCTX = [tUpCTX; (cell2mat(allCTX{recNum,s}.evtPeak) - cell2mat(allCTX{recNum,s}.evtBegin))*dt];
            tDownCTX = [tDownCTX; (cell2mat(allCTX{recNum,s}.evtEnd) - cell2mat(allCTX{recNum,s}.evtPeak))*dt];
            evtFreqCTX = [evtFreqCTX; allCTX{recNum,s}.meanEvtFreq];
            evtAmpliCTX = [evtAmpliCTX; cell2mat(allCTX{recNum,s}.evtAmpli)];
            actNumCTX = [actNumCTX; numel(allCTX{recNum,s}.allActiveAstros)];
            totNumCTX = [totNumCTX; allCTX{recNum,s}.nAstros];
            % *** Remove null events *** %
            evtdurCTX = evtdurCTX(evtdurCTX~=0);
            tUpCTX = tUpCTX(evtdurCTX~=0);
            tDownCTX = tDownCTX(evtdurCTX~=0);
            evtAmpliCTX = evtAmpliCTX(evtdurCTX~=0);
            evtFreqCTX = evtFreqCTX(evtFreqCTX~=0);
            interEvtCTX = interEvtCTX(interEvtCTX~=0);
            % Sub events (ctrl CTX)
            numSub = [];
            for i = 1:length(allCTX{recNum,s}.evtSubPeaks)
                numSub = [numSub; (cell2mat(allCTX{recNum,s}.evtSubPeaks(i)))];
            end
            numSubEvtCTX = [numSubEvtCTX; numSub];
            freqSubEvtCTX = [freqSubEvtCTX; cell2mat(allCTX{recNum,s}.evtSubPeaksFreq)];
        end
    end
end

% Concatenate values for all sessions CA1 ----------------------------- %
% Initiate (CA1)
tUpCA1 = [];
tDownCA1 = [];
evtdurCA1 = [];
evtFreqCA1 = [];
evtAmpliCA1 = [];
numSubEvtCA1 = [];
freqSubEvtCA1 = [];
interEvtCA1 = [];
actNumCA1 = [];
totNumCA1 = [];

% concatenate (CA1)
for recNum = recList
    dt = 0.5;
    for s = 1:length(allCTX(recNum,:))
        %%%% ALL CTRL SESSIONS CA1 %%%%
        if ~isempty(allCA1{recNum,s})
            evtdurCA1 = [evtdurCA1; (cell2mat(allCA1{recNum,s}.evtEnd) - cell2mat(allCA1{recNum,s}.evtBegin))*dt];
            for e = 1:length(allCA1{recNum,s}.evtBegin)
                if length(cell2mat(allCA1{recNum,s}.evtBegin(e)))>=2
                    beg = cell2mat(allCA1{recNum,s}.evtBegin(e));
                    nd = cell2mat(allCA1{recNum,s}.evtEnd(e));
                    interEvtCA1 = [interEvtCA1; (beg(2:end)-nd(1:end-1))*dt];
                end
            end
            tUpCA1 = [tUpCA1; (cell2mat(allCA1{recNum,s}.evtPeak) - cell2mat(allCA1{recNum,s}.evtBegin))*dt];
            tDownCA1 = [tDownCA1; (cell2mat(allCA1{recNum,s}.evtEnd) - cell2mat(allCA1{recNum,s}.evtPeak))*dt];
            evtFreqCA1 = [evtFreqCA1; allCA1{recNum,s}.meanEvtFreq];
            evtAmpliCA1 = [evtAmpliCA1; cell2mat(allCA1{recNum,s}.evtAmpli)];
            actNumCA1 = [actNumCA1; numel(allCA1{recNum,s}.allActiveAstros)];
            totNumCA1 = [totNumCA1; allCA1{recNum,s}.nAstros];
            % *** Remove null events *** %
            evtdurCA1 = evtdurCA1(evtdurCA1~=0);
            tUpCA1 = tUpCA1(evtdurCA1~=0);
            tDownCA1 = tDownCA1(evtdurCA1~=0);
            evtAmpliCA1 = evtAmpliCA1(evtdurCA1~=0);
            evtFreqCA1 = evtFreqCA1(evtFreqCA1~=0);
            interEvtCA1 = interEvtCA1(interEvtCA1~=0);
            % Sub events (CA1 ctrl)
            numSub = [];
            for i = 1:length(allCA1{recNum,s}.evtSubPeaks)
                numSub = [numSub; (cell2mat(allCA1{recNum,s}.evtSubPeaks(i)))];
            end
            numSubEvtCA1 = [numSubEvtCA1; numSub];
            freqSubEvtCA1 = [freqSubEvtCA1; cell2mat(allCA1{recNum,s}.evtSubPeaksFreq)];
        end
    end
end

%% Number of peaks per event
cmap = hsv;
CTXcol = 10;
CA1col = 50;
% Max number of peaks over all sessions
maxPeak = max([max(numSubEvtCTX), max(numSubEvtCA1)]);
% ------------------------ Concatenate for bar plot --------------------- %
barPeaks = zeros(maxPeak,2);
for i = 1:maxPeak
    barPeaks(i,:) = [numel(numSubEvtCTX(numSubEvtCTX == i))/numel(numSubEvtCTX) ...
        numel(numSubEvtCA1(numSubEvtCA1 == i))/numel(numSubEvtCA1) ...
        ];
end
figure
b = bar(barPeaks(1:5,:),'FaceColor','flat','EdgeColor','none');
b(1).FaceColor = cmap(CTXcol,:);
b(2).FaceColor = cmap(CA1col,:);
title('Number of peaks within events')
xlabel('Number of peaks')
ylabel('Probability')
legend(['CTX: ' num2str(mean(numSubEvtCTX),'%4.2f') '\pm' ...
    num2str(std(numSubEvtCTX)/sqrt(numel(numSubEvtCTX)),'%4.2f')],...
    ['CA1: ' num2str(mean(numSubEvtCA1),'%4.2f') '\pm'...
    num2str(std(numSubEvtCA1)/sqrt(numel(numSubEvtCA1)),'%4.2f')])
% Stats ks2
[h_ctx_CA1,p_ctx_CA1] = kstest2(numSubEvtCTX,numSubEvtCA1);

%% Events durations -- Histograms
figure
% CTX
histogram(evtdurCTX,'FaceColor',cmap(CTXcol,:),'DisplayStyle','bar','EdgeColor','none',...
    'binEdges',binEdg,'normalization','probability','FaceAlpha',0.75)
hold on 
% CA1
histogram(evtdurCA1,'FaceColor','none','DisplayStyle','stairs',...
    'EdgeColor',cmap(CA1col,:),...
    'binEdges',binEdg,'normalization','probability','LineWidth',1.5)

legend('CTX','CA1')
xlabel('Duration (s)')
ylabel('Probability')
title('Event durations')

% Durations -- Violins
maxLengthDur = max([length(evtdurCTX),length(evtdurCA1)]);

durations = [[evtdurCTX; nan(maxLengthDur-length(evtdurCTX),1)]...
    [evtdurCA1; nan(maxLengthDur-length(evtdurCA1),1)]...
    ];
figure
violinsDur = violinplot(durations,{'CTX','CA1'},...
    'ShowMean', true, 'ShowData', false, 'ViolinAlpha', 1);
violinsDur(1).ViolinColor = cmap(CTXcol,:);
violinsDur(2).ViolinColor = cmap(CA1col,:);
title('Event durations')
ylabel('Duration (s)')
ylim([0,30])
% Add mean values on the plot
text(1.05,mean(evtdurCTX),[num2str(mean(evtdurCTX),'%4.2f') '\pm' num2str(std(evtdurCTX),'%4.2f') 's'])
text(2.05,mean(evtdurCA1),[num2str(mean(evtdurCA1),'%4.2f') '\pm' num2str(std(evtdurCA1),'%4.2f') 's'])
%% Amplitude -- Histograms
binEdg = 0:0.1:15;
figure
% CTX
histogram(evtAmpliCTX,'FaceColor',cmap(CTXcol,:),'DisplayStyle','bar','EdgeColor','none',...
    'binEdges',binEdg,'normalization','probability','FaceAlpha',0.75)
hold on 
% CA1
histogram(evtAmpliCA1,'EdgeColor',cmap(CA1col,:),'DisplayStyle','stairs',...
    'FaceColor','none',...
    'binEdges',binEdg,'normalization','probability','LineWidth',1.5)

legend('CTX','CA1')
xlabel('Amplitude')
ylabel('Probability')
title('Events amplitudes')

% Amplitudes -- Violins
maxLengthAmpli = max([length(evtAmpliCTX),length(evtAmpliCA1)]);
amplitudes = [[evtAmpliCTX; nan(maxLengthAmpli-length(evtAmpliCTX),1)]...
    [evtAmpliCA1; nan(maxLengthAmpli-length(evtAmpliCA1),1)]...
    ];
figure
violinsAmpli = violinplot(amplitudes,{'VLPO F','VLPO M'},...
    'ShowMean',true, 'ShowData', false, 'ViolinAlpha', 1);
violinsAmpli(1).ViolinColor = cmap(CTXcol,:);
violinsAmpli(2).ViolinColor = cmap(CA1col,:);
title('Event amplitudes')
ylabel('Amplitudes')
ylim([0,15])
% Add mean values on the plot
text(1.05,mean(evtAmpliCTX),[num2str(mean(evtAmpliCTX),'%4.2f') '\pm' num2str(std(evtAmpliCTX),'%4.2f')])
text(2.05,mean(evtAmpliCA1),[num2str(mean(evtAmpliCA1),'%4.2f') '\pm' num2str(std(evtAmpliCA1),'%4.2f')])
%% Frequency
binEdg = 0:0.025:0.5;
figure
% CTX
histogram(evtFreqCTX,'FaceColor',cmap(CTXcol,:),'DisplayStyle','bar','EdgeColor','none',...
    'binEdges',binEdg,'normalization','probability','FaceAlpha',0.75)
hold on
% CA1
histogram(evtFreqCA1,'EdgeColor',cmap(CA1col,:),...
    'DisplayStyle','stairs','FaceColor','none',...
    'binEdges',binEdg,'normalization','probability','LineWidth',1.5)
legend('CTX','CA1')
xlabel('Frequency (Hz)')
ylabel('Probability')

% Frequency -- Violins
maxLengthFreq = max([length(evtFreqCTX),length(evtFreqCA1)]);

frequences = [[evtFreqCTX; nan(maxLengthFreq-length(evtFreqCTX),1)]...
    [evtFreqCA1; nan(maxLengthFreq-length(evtFreqCA1),1)]...
    ];
figure
violinsFreq = violinplot(frequences,{'VLPO F','VLPO M'},...
    'ShowMean',true, 'ShowData', false, 'ViolinAlpha', 1);
violinsFreq(1).ViolinColor = cmap(CTXcol,:);
violinsFreq(2).ViolinColor = cmap(CA1col,:);
title('Event frequencies')
ylabel('Frequency (Hz)')
ylim([0,0.7])
% Add mean values on the plot
text(1.05,mean(evtFreqCTX),[num2str(mean(evtFreqCTX),'%4.2f') '\pm' num2str(std(evtFreqCTX),'%4.2f') 'Hz'])
text(2.05,mean(evtFreqCA1),[num2str(mean(evtFreqCA1),'%4.2f') '\pm' num2str(std(evtFreqCA1),'%4.2f') 'Hz'])
%% Correlation curves
corLevelMatCTX = [];
for recNum = 1:2:3
    for s = 1:5
        %%%% CONCATENATE ALL SESSIONS CTX %%%%
        if ~isempty(allCTX{recNum,s})
            corLevelMatCTX = [corLevelMatCTX; allCTX{recNum,s}.corrCurve];
        end
    end
end

corLevelMatCA1 = [];
for recNum = 1:2:3
    for s = 1:5
        %%%% CONCATENATE ALL SESSIONS CA1 %%%%
        if ~isempty(allCA1{recNum,s}) 
            corLevelMatCA1 = [corLevelMatCA1; allCA1{recNum,s}.corrCurve];
        end
    end
end

% Correlation plots mean \pm std
corThresh = 0:0.01:1;
figure
% CTX
plot(corThresh,mean(corLevelMatCTX,1),'color',cmap(CTXcol,:),'LineWidth',1.5)
hold on
patch([corThresh fliplr(corThresh)], ...
    [mean(corLevelMatCTX,1)+std(corLevelMatCTX,1)/sqrt(size(corLevelMatCTX,1))...
    fliplr(mean(corLevelMatCTX,1)-std(corLevelMatCTX,1)/sqrt(size(corLevelMatCTX,1)))],...
    cmap(CTXcol,:),'EdgeColor','none','FaceAlpha',0.3,'HandleVisibility','off')
areaCTX = trapz(corThresh,mean(corLevelMatCTX,1));
areaSupCTX = trapz(corThresh,mean(corLevelMatCTX,1)+...
    std(corLevelMatCTX,1)/sqrt(size(corLevelMatCTX,1)))-areaCTX;
areaMinCTX = areaCTX-trapz(corThresh,mean(corLevelMatCTX,1)-...
    std(corLevelMatCTX,1)/sqrt(size(corLevelMatCTX,1)));
% CA1
plot(corThresh,mean(corLevelMatCA1,1),'color',cmap(CA1col,:),'LineWidth',1.5)
hold on
patch([corThresh fliplr(corThresh)], ...
    [mean(corLevelMatCA1,1)+std(corLevelMatCA1,1)/sqrt(size(corLevelMatCA1,1))...
    fliplr(mean(corLevelMatCA1,1)-std(corLevelMatCA1,1)/sqrt(size(corLevelMatCA1,1)))],...
    cmap(CA1col,:),'EdgeColor','none','FaceAlpha',0.3,'HandleVisibility','off')
areaCA1 = trapz(corThresh,mean(corLevelMatCA1,1));
areaSupCA1 = trapz(corThresh,mean(corLevelMatCA1,1)+...
    std(corLevelMatCA1,1)/sqrt(size(corLevelMatCA1,1)))-areaCA1;
areaMinCA1 = areaCA1-trapz(corThresh,mean(corLevelMatCA1,1)-...
    std(corLevelMatCA1,1)/sqrt(size(corLevelMatCA1,1)));
legend('CTX','CA1')
ylim([0,1])
xlabel('correlation threshold')
ylabel('Proportion of astrocytes (%)')
title('Proportion of astrocytes over correlation threshold')

[h_ctx_CA1,p_ctx_CA1] = kstest2(corLevelMatCTX(:),corLevelMatCA1(:));