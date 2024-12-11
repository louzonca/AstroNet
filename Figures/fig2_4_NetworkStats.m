% This script concatenates all the processed data contained in the
% allCTX and allCA1 mat files in Matlab arrays to make the plots of the
% network statistics of the events (Figs. 2 and 4 in the manuscript)
% --- Load data --------------------------------------------------------- %
load('Data\allCTX.mat')
load('Data\allCA1.mat')

%% Global event data
recList = 1:2:5;
segmentations = 1:1:20;

% Initiate CTX
pathLengthsCTX = [];
numActiveAstroPerEventCTX = [];
activeAreaCTX = [];
normPathsCTX = [];
normActiveAreaCTX = [];
meanDegCTX = nan(5,8,length(segmentations));
varDegCTX = nan(5,8,length(segmentations));
numHighConnNodesCTX = nan(5,8,length(segmentations));
numHighConnNodesNormCTX = nan(5,8,length(segmentations));

% concatenate (CTX)
for expType = recList
    dt = 0.5;
    for s = 1:5
        %%%% ALL CTRL SESSIONS CTX %%%%
        if ~isempty(allCTX{expType,s})
            % Raw global event values CTX
            pathLengthsCTX = [pathLengthsCTX; allCTX{expType,s}.pathLengths];
            activeAreaCTX = [activeAreaCTX; allCTX{expType,s}.activeAreas];            
            numActiveAstroPerEventCTX = [numActiveAstroPerEventCTX; allCTX{expType,s}.activeAstrosEvt];
            % Normalized global event values CTX
            normPathsCTX = [normPathsCTX; allCTX{expType,s}.pathLengths./allCTX{expType,s}.activeAstrosEvt];
            normActiveAreaCTX = [normActiveAreaCTX; allCTX{expType,s}.activeAreas./allCTX{expType,s}.activeAstrosEvt]; 
            % Remove 0s CTX
            pathLengthsCTX = pathLengthsCTX(pathLengthsCTX~=0);
            numActiveAstroPerEventCTX = numActiveAstroPerEventCTX(numActiveAstroPerEventCTX~=0);
            activeAreaCTX = activeAreaCTX(activeAreaCTX~=0 & ~isnan(activeAreaCTX));
            % Graphs stats CTX
            meanDegCTX(expType,s,:) = allCTX{expType,s}.meanNodeDeg;
            varDegCTX(expType,s,:) = allCTX{expType,s}.varNodeDeg;
            numHighConnNodesCTX(expType,s,:) = allCTX{expType,s}.numHighConnNodes;
            numHighConnNodesNormCTX(expType,s,:) = allCTX{expType,s}.numHighConnNodes/length(allCTX{expType,s}.allActiveAstros);
              %  allCTX{expType,s}.propHighConnNodes;
        end
    end
end

% CA1
% Initiate CA1
pathLengthsCA1 = [];
normPathsCA1 = [];
numActiveAstroPerEventCA1 = [];
activeAreaCA1 = [];
normActiveAreaCA1 = [];
meanDegCA1 = nan(5,12,length(segmentations));
varDegCA1 = nan(5,12,length(segmentations));
numHighConnNodesCA1 = nan(5,12,length(segmentations));
numHighConnNodesNormCA1 = nan(5,12,length(segmentations));

% concatenate (CA1)
for expType = recList
    dt = 0.5;
    for s = 1:5
        %%%% ALL CTRL SESSIONS CA1 %%%%
        if ~isempty(allCA1{expType,s})
            % Raw global event values CA1
            pathLengthsCA1 = [pathLengthsCA1; allCA1{expType,s}.pathLengths];
            activeAreaCA1 = [activeAreaCA1; allCA1{expType,s}.activeAreas];            
            numActiveAstroPerEventCA1 = [numActiveAstroPerEventCA1; allCA1{expType,s}.activeAstrosEvt];
            % Normalized global event values CA1
            normPathsCA1 = [normPathsCA1; allCA1{expType,s}.pathLengths./allCA1{expType,s}.activeAstrosEvt];
            normActiveAreaCA1 = [normActiveAreaCA1; allCA1{expType,s}.activeAreas./allCA1{expType,s}.activeAstrosEvt]; 
            % Remove 0s CA1
            pathLengthsCA1 = pathLengthsCA1(pathLengthsCA1~=0);
            numActiveAstroPerEventCA1 = numActiveAstroPerEventCA1(numActiveAstroPerEventCA1~=0);
            activeAreaCA1 = activeAreaCA1(activeAreaCA1~=0 & ~isnan(activeAreaCA1));
            % Graphs stats CA1
            meanDegCA1(expType,s,:) = allCA1{expType,s}.meanNodeDeg;
            varDegCA1(expType,s,:) = allCA1{expType,s}.varNodeDeg;
            numHighConnNodesCA1(expType,s,:) = allCA1{expType,s}.numHighConnNodes;
            numHighConnNodesNormCA1(expType,s,:) = allCA1{expType,s}.numHighConnNodes/length(allCA1{expType,s}.allActiveAstros);              
%allCA1{expType,s}.propHighConnNodes;
        end
    end
end

%% Path lengths -- Violins
cmap = hsv;
CTXcol = 10;
CA1col = 50;

maxLengthPaths = max([length(pathLengthsCTX),length(pathLengthsCA1)]);

paths = [[pathLengthsCTX; nan(maxLengthPaths-length(pathLengthsCTX),1)]...
    [pathLengthsCA1; nan(maxLengthPaths-length(pathLengthsCA1),1)]...
    ];
figure
violinsPath = violinplot(paths,{'CTX', 'CA1'}...
    ,'ShowMean',true, 'showData', false, 'ViolinAlpha', 1);
violinsPath(1).ViolinColor = cmap(CTXcol,:);
violinsPath(2).ViolinColor = cmap(CA1col,:);
title('Path length')
ylabel('Length (um)')
ylim([0,20000])
% Add mean values on the plot
text(1.05,mean(pathLengthsCTX),[num2str(mean(pathLengthsCTX),'%4.2f') '\pm' num2str(std(pathLengthsCTX),'%4.2f') 'um'])
text(2.05,mean(pathLengthsCA1),[num2str(mean(pathLengthsCA1),'%4.2f') '\pm' num2str(std(pathLengthsCA1),'%4.2f') 'um'])
% ks2 stats

%% Path length/number of active astros
% Path lengths
maxLengthNormedPaths = max([length(normPathsCTX),length(normPathsCA1)]);

normedPaths = [[normPathsCTX; nan(maxLengthNormedPaths-length(normPathsCTX),1)]...
    [normPathsCA1; nan(maxLengthNormedPaths-length(normPathsCA1),1)]...
    ];
figure
violinsPath = violinplot(normedPaths,{'CTX', 'CA1'},...
    'ShowMean',true, 'showData', false, 'ViolinAlpha', 1);
violinsPath(1).ViolinColor = cmap(CTXcol,:);
violinsPath(2).ViolinColor = cmap(CA1col,:);
title('Path length/number of active astrocytes')
ylabel('Length (um)')
ylim([200,600])
% Add mean values on the plot
normPathsCTX = normPathsCTX(~isnan(normPathsCTX));
normPathsCA1 = normPathsCA1(~isnan(normPathsCA1));
text(1.05,mean(normPathsCTX),[num2str(mean(normPathsCTX),'%4.2f') '\pm' num2str(std(normPathsCTX),'%4.2f') 'um'])
text(2.05,mean(normPathsCA1),[num2str(mean(normPathsCA1),'%4.2f') '\pm' num2str(std(normPathsCA1),'%4.2f') 'um'])

% ks2 stats
[h_ctx_CA1,p_ctx_CA1] = kstest2(normPathsCTX,normPathsCA1);

%% Active areas -- Violins
maxLengthAreas = max([length(activeAreaCTX),length(activeAreaCA1)]);

areas = [[activeAreaCTX; nan(maxLengthAreas-length(activeAreaCTX),1)]...
    [activeAreaCA1; nan(maxLengthAreas-length(activeAreaCA1),1)]...
    ];
figure
violinsArea = violinplot(areas,{'CTX', 'CA1'},...
    'ShowMean',true, 'showData', false, 'ViolinAlpha', 1);
violinsArea(1).ViolinColor = cmap(CTXcol,:);
violinsArea(2).ViolinColor = cmap(CA1col,:);
title('Active area')
ylabel('Area (um2)')
% Add mean values on the plot
text(1.05,mean(activeAreaCTX),[num2str(mean(activeAreaCTX),'%4.2f') '\pm' num2str(std(activeAreaCTX),'%4.2f') 'um2'])
text(2.05,mean(activeAreaCA1),[num2str(mean(activeAreaCA1),'%4.2f') '\pm' num2str(std(activeAreaCA1),'%4.2f') 'um2'])

% ks2 stats
[h_ctx_CA1,p_ctx_CA1] = kstest2(activeAreaCTX,activeAreaCA1);

%% Active areas/number of active astrocytes -- Violins
maxLengthNormAreas = max([length(normActiveAreaCTX),length(normActiveAreaCA1)]);

normAreas = [[normActiveAreaCTX; nan(maxLengthNormAreas-length(normActiveAreaCTX),1)]...
    [normActiveAreaCA1; nan(maxLengthNormAreas-length(normActiveAreaCA1),1)]...
    ];
figure
violinsArea = violinplot(normAreas,{'CTX', 'CA1','VLPO'},...
    'ShowMean',true, 'showData', false, 'ViolinAlpha', 1);
violinsArea(1).ViolinColor = cmap(CTXcol,:);
violinsArea(2).ViolinColor = cmap(CA1col,:);
title('Active area/number of active astrocytes')
ylabel('Area (um2)')
ylim([6000,60000])
% Add mean values on the plot
normActiveAreaCA1 = normActiveAreaCA1(~isnan(normActiveAreaCA1));
normActiveAreaCTX = normActiveAreaCTX(~isnan(normActiveAreaCTX));
text(2.05,mean(normActiveAreaCA1),[num2str(mean(normActiveAreaCA1),'%4.2f') '\pm' num2str(std(normActiveAreaCA1),'%4.2f') 'um2'])
text(1.05,mean(normActiveAreaCTX),[num2str(mean(normActiveAreaCTX),'%4.2f') '\pm' num2str(std(normActiveAreaCTX),'%4.2f') 'um2'])

% ks2 stats
[h_ctx_CA1,p_ctx_CA1] = kstest2(normActiveAreaCTX,normActiveAreaCA1);

%% Number of active cells -- Violins
maxLengthActives = max([length(numActiveAstroPerEventCA1),...
    length(numActiveAstroPerEventCTX)]);

actives = [[numActiveAstroPerEventCTX; nan(maxLengthActives-length(numActiveAstroPerEventCTX),1)]...
    [numActiveAstroPerEventCA1; nan(maxLengthActives-length(numActiveAstroPerEventCA1),1)]...
    ];
figure
violinsActives = violinplot(actives,{'CTX', 'CA1'},...
    'ShowMean',true, 'showData', false, 'ViolinAlpha', 1);
violinsActives(1).ViolinColor = cmap(CTXcol,:);
violinsActives(2).ViolinColor = cmap(CA1col,:);
title('Active astrocytes per global event')
ylabel('Number')
% Add mean values on the plot
text(2.05,mean(numActiveAstroPerEventCA1),[num2str(mean(numActiveAstroPerEventCA1),'%4.2f') '\pm' num2str(std(numActiveAstroPerEventCA1),'%4.2f')])
text(1.05,mean(numActiveAstroPerEventCTX),[num2str(mean(numActiveAstroPerEventCTX),'%4.2f') '\pm' num2str(std(numActiveAstroPerEventCTX),'%4.2f')])

% ks2 stats
[h_ctx_CA1,p_ctx_CA1] = kstest2(numActiveAstroPerEventCTX,numActiveAstroPerEventCA1);

%% Number of highly connected nodes
maxLengthHC = max([length(numHighConnNodesCA1(:)),...
    length(numHighConnNodesCTX(:))]);

HC = [[numHighConnNodesCTX(:); nan(maxLengthHC-length(numHighConnNodesCTX(:)),1)]...
    [numHighConnNodesCA1(:); nan(maxLengthHC-length(numHighConnNodesCA1(:)),1)]...
    ];
figure
violinsHC = violinplot(HC,{'CTX', 'CA1'},...
    'ShowMean',true, 'showData', false, 'ViolinAlpha', 1);
violinsHC(1).ViolinColor = cmap(CTXcol,:);
violinsHC(2).ViolinColor = cmap(CA1col,:);
title('Highy connected astrocytes')
ylabel('Number')
hold on
% Add mean values on the plot
numHighConnNodesCA1 = numHighConnNodesCA1(~isnan(numHighConnNodesCA1));
numHighConnNodesCTX = numHighConnNodesCTX(~isnan(numHighConnNodesCTX));
text(2.05,mean(numHighConnNodesCA1(:)),[num2str(mean(numHighConnNodesCA1(:)),'%4.2f') '\pm'...
    num2str(std(numHighConnNodesCA1(:)),'%4.2f')])
text(1.05,mean(numHighConnNodesCTX(:)),[num2str(mean(numHighConnNodesCTX(:)),'%4.2f') '\pm'...
    num2str(std(numHighConnNodesCTX(:)),'%4.2f')])

% ks2 stats
[h_ctx_CA1,p_ctx_CA1] = kstest2(numHighConnNodesCTX,numHighConnNodesCA1);

%% Number of highly connected nodes/totalnumber of active astros
maxLengthHCnorm = max([length(numHighConnNodesNormCA1(:)),...
    length(numHighConnNodesNormCTX(:))]);

HCnorm = [[numHighConnNodesNormCTX(:); nan(maxLengthHCnorm-length(numHighConnNodesNormCTX(:)),1)]...
    [numHighConnNodesNormCA1(:); nan(maxLengthHCnorm-length(numHighConnNodesNormCA1(:)),1)]...
    ];
figure
violinsHCnorm = violinplot(HCnorm,{'CTX', 'CA1'},...
    'ShowMean',true, 'showData', false, 'ViolinAlpha', 1);
violinsHCnorm(1).ViolinColor = cmap(CTXcol,:);
violinsHCnorm(2).ViolinColor = cmap(CA1col,:);
title('Highy connected astrocytes (normalized)')
ylabel('Number')
% Add mean values on the plot
numHighConnNodesNormCA1 = numHighConnNodesNormCA1(~isnan(numHighConnNodesNormCA1));
numHighConnNodesNormCTX = numHighConnNodesNormCTX(~isnan(numHighConnNodesNormCTX));
text(2.05,mean(numHighConnNodesNormCA1(:)),[num2str(mean(numHighConnNodesNormCA1(:)),'%4.2f') '\pm'...
    num2str(std(numHighConnNodesNormCA1(:)),'%4.2f')])
text(1.05,mean(numHighConnNodesNormCTX(:)),[num2str(mean(numHighConnNodesNormCTX(:)),'%4.2f') '\pm'...
    num2str(std(numHighConnNodesNormCTX(:)),'%4.2f')])

% ks2 stats
[h_ctx_CA1,p_ctx_CA1] = kstest2(numHighConnNodesNormCTX,numHighConnNodesNormCA1);

%% Mean node degree
maxLengthNodeDeg = max([length(meanDegCA1(:)),length(meanDegCTX(:))]);

nodeDeg = [[meanDegCTX(:); nan(maxLengthNodeDeg-length(meanDegCTX(:)),1)]...
    [meanDegCA1(:); nan(maxLengthNodeDeg-length(meanDegCA1(:)),1)]...
    ];
figure
violinsNodeDeg = violinplot(nodeDeg,{'CTX', 'CA1'},...
    'ShowMean',true, 'showData', false, 'ViolinAlpha', 1);
violinsNodeDeg(1).ViolinColor = cmap(CTXcol,:);
violinsNodeDeg(2).ViolinColor = cmap(CA1col,:);
title('Mean astrocytes connection degree')
ylabel('Number')
% Add mean values
meanDegCA1 = meanDegCA1(~isnan(meanDegCA1));
meanDegCTX = meanDegCTX(~isnan(meanDegCTX));
text(2.05,mean(meanDegCA1(:)),[num2str(mean(meanDegCA1(:)),'%4.2f') '\pm' num2str(std(meanDegCA1(:)),'%4.2f')])
text(1.05,mean(meanDegCTX(:)),[num2str(mean(meanDegCTX(:)),'%4.2f') '\pm' num2str(std(meanDegCTX(:)),'%4.2f')])

% ks2 stats
[h_ctx_CA1,p_ctx_CA1] = kstest2(meanDegCTX,meanDegCA1);

%% Node degree sorted distribution
meanDegCA1sorted = sort(meanDegCA1);
meanDegCTXsorted = sort(meanDegCTX);
figure
hold on
histogram(meanDegCTXsorted,8,'DisplayStyle','bar','FaceColor',cmap(CTXcol,:),'EdgeColor','none','FaceAlpha',1,'normalization','probability')
histogram(meanDegCA1sorted,8,'DisplayStyle','stairs','EdgeColor',cmap(CA1col,:),'LineWidth',2,'FaceColor','none','normalization','probability')
legend('CTX','CA1')

%% Concatenate data for paths steps plots
% CTX %
recList = 1:2:5;
pathStepsCTX = [];
pathStepsTimeCTX = [];

for expType = recList
    dt = 0.5;
    for s = 1:5
        if ~isempty(allCTX{expType,s})
            for pl=1:length(allCTX{expType,s}.pathSteps)
                pathStepsCTX = [pathStepsCTX; allCTX{expType,s}.pathSteps{pl}];
                pathStepsTimeCTX = [pathStepsTimeCTX; allCTX{expType,s}.pathStepsTimes{pl}*dt];
            end
        end
    end
end

% CA1 %
recList = 1:2:5;
pathStepsCA1 = [];
pathStepsTimeCA1 = [];

for expType = recList
     for s = 1:5
         if ~isempty(allCA1{expType,s})
            for pl=1:length(allCA1{expType,s}.pathSteps)
                pathStepsCA1 = [pathStepsCA1; allCA1{expType,s}.pathSteps{pl}];
                pathStepsTimeCA1 = [pathStepsTimeCA1; allCA1{expType,s}.pathStepsTimes{pl}*dt];
            end
         end
    end
end

%% Paths steps -- Histograms
cmap = hsv;
CTXcol = 10;
CA1col = 50;
figure
hold on
histogram(pathStepsCTX,30,'DisplayStyle','bar','FaceColor',cmap(CTXcol,:),'EdgeColor','none','FaceAlpha',1,'normalization','probability')
histogram(pathStepsCA1,30,'DisplayStyle','stairs','EdgeColor',cmap(CA1col,:),'LineWidth',2,'FaceColor','none','normalization','probability')
legend('CTX','CA1')
xlabel('Steps length')
ylabel('Probability')
% ks2 stats
[h_ctx_CA1,p_ctx_CA1] = kstest2(pathStepsCTX,pathStepsCA1);

% Mean \pm sd
mCTX = mean(pathStepsCTX);
sCTX = std(pathStepsCTX)/sqrt(length(pathStepsCTX));

mCA1 = mean(pathStepsCA1);
sCA1 = std(pathStepsCA1)/sqrt(length(pathStepsCA1));

%% Paths steps durations -- Histograms
bins = 0:max(pathStepsTimeCTX)/46:max(pathStepsTimeCTX);
figure
hold on
hCTX = histogram(abs(pathStepsTimeCTX),40,'DisplayStyle','bar','FaceColor',...
    cmap(CTXcol,:),'EdgeColor','none','FaceAlpha',0.6,'normalization','probability')
ctxVal = hCTX.Values;
ctxBin = hCTX.BinEdges(2:end);
a = 0.1381; % fit was done with cftool
b = -1.176;
plot(ctxBin,a.*ctxBin.^b,'--','color',cmap(CTXcol,:),'LineWidth',1.5)
l1 = [num2str(a) 'x^' num2str(b)];
hCA1 = histogram(abs(pathStepsTimeCA1),40,'DisplayStyle','bar','FaceColor',...
    cmap(CA1col,:),'EdgeColor','none','FaceAlpha',0.5,'normalization','probability')
CA1Val = hCA1.Values;
CA1Bin = hCA1.BinEdges(2:end);
a = 0.09937;
b = -1.779;
l2 = [num2str(a) 'x^' num2str(b)];
plot(CA1Bin,a.*CA1Bin.^b,'--','color',cmap(CA1col-2,:),'LineWidth',1.5)

legend('CTX',l1,'CA1',l2)
xlabel('Steps duration (s)')
ylabel('Probability')
xlim([0,5])
ylim([0,0.5])
% ks2 stats
[h_ctx_CA1,p_ctx_CA1] = kstest2(pathStepsTimeCTX,pathStepsTimeCA1);

%% Distribution of edges weights
edgesWeightsCTX = [];
for expType = recList
    for s = 1:5
        if ~isempty(allCTX{expType,s})
            edgesWeightsCTX = [edgesWeightsCTX; allCTX{expType,s}.conMat(:)];
        end
    end
end

region = 2;
edgesWeightsCA1 = [];
for expType = recList
    for s = 1:5
        if ~isempty(allCA1{expType,s})
            edgesWeightsCA1 = [edgesWeightsCA1; allCA1{expType,s}.conMat(:)];
        end
    end
end

%% Edges weights -- Violins 
maxLengthEW = max([length(edgesWeightsCTX),length(edgesWeightsCA1)]);

ew = [[edgesWeightsCTX; nan(maxLengthEW-length(edgesWeightsCTX),1)]...
    [edgesWeightsCA1; nan(maxLengthEW-length(edgesWeightsCA1),1)]...
    ];
figure
violinsPath = violinplot(ew,{'CTX', 'CA1'}...
    ,'ShowMean',true, 'showData', false, 'ViolinAlpha', 1);
violinsPath(1).ViolinColor = cmap(CTXcol,:);
violinsPath(2).ViolinColor = cmap(CA1col,:);
title('Edges weights')
ylabel('weigth')
%ylim([0,1])
% Add mean values on the plot
text(1.05,mean(edgesWeightsCTX),[num2str(mean(edgesWeightsCTX),'%4.2f') '\pm' num2str(std(edgesWeightsCTX),'%4.2f') 'um'])
text(2.05,mean(edgesWeightsCA1),[num2str(mean(edgesWeightsCA1),'%4.2f') '\pm' num2str(std(edgesWeightsCA1),'%4.2f') 'um'])
% ks2 stats
[h_ctx_CA1,p_ctx_CA1] = kstest2(edgesWeightsCTX,edgesWeightsCA1);

%% Edges weights -- Histograms
bins = 1:2:20;
figure
hold on
histogram(edgesWeightsCTX(edgesWeightsCTX>0),'BinEdges',bins,'DisplayStyle','bar',...
    'FaceColor',cmap(CTXcol,:),'EdgeColor','none','FaceAlpha',1,'normalization','probability')
histogram(edgesWeightsCA1(edgesWeightsCA1>0),'BinEdges',bins,'DisplayStyle','stairs',...
    'EdgeColor',cmap(CA1col,:),'LineWidth',2,'FaceColor','none','normalization','probability')
legend('CTX','CA1')
xlabel('Edges weights')
ylabel('Probability')
