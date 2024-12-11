% Load data
path = 'C:\Users\louzo\Dropbox\IBENS\05_Giampaolo_calciumNetworks\01_code\data\timelapse\Vidéos Imagerie calcique aout 2024\';

dates = {' 06 ', ' 07 ', ' 08 ',' 09 '};
regions = {'X1826 VLPO M1 F', 'X1827 VLPO M1 F', 'X1828 VLPO CA1 F', 'X1829 VLPO CA1 F'};
nums = {[1], [7], [3, 5], [5, 10, 10.2, 12]};
regName = {' Cx Mtr ', ' Cx Mtr ', ' CA1 ', ' CA1 '};

ssimCTX = [];
aucCTX = [];
ksCTX = [];
ssimCA1 = [];
aucCA1 = [];
ksCA1 = [];

propHCCTX = [];
numHCCTX = [];
meanDegCTX = [];
varDegCTX = [];

propHCCA1 = [];
numHCCA1 = [];
meanDegCA1 = [];
varDegCA1 = [];

areaCTX = [];
numActivCTX = [];
evPathLengthCTX = [];

areaCA1 = [];
numActivCA1 = [];
evPathLengthCA1 = [];

for n = 1:4
    date = dates{n};
    region = regions{n};
    numbers = nums{n};
    reg = regName{n};
    
    for num = numbers
        disp('Detecting astros')
        numpt = 1;
        v_1 = VideoReader([path '2024 08' date region '\2024 08' date 'F' reg num2str(num) ' pt ' ...
            num2str(numpt) '.wmv']); % Cx Mtr %CA1
        data1 = read(v_1,[1 Inf]);
        data1 = squeeze(data1(:,:,1,:));

        numpt = 2;
        v_2 = VideoReader([path '2024 08' date region '\2024 08' date 'F' reg num2str(num) ' pt ' ...
            num2str(numpt) '.wmv']);
        data2 = read(v_2,[1 Inf]);
        data2 = squeeze(data2(:,:,1,:));

        %implay(data1)

        % Detect astrocytes
        sumImg1 = sum(data1(:,:,:),3);
        sumImg2 = sum(data2(:,:,:),3);
        figure 
        subplot(1,2,1)
        colormap 
        imshow(sumImg1*2,[min(sumImg1(:)),max(sumImg1(:))])
        subplot(1,2,2)
        colormap 
        imshow(sumImg2*2,[min(sumImg2(:)),max(sumImg2(:))])
        suptitle(['Summed image - session ' region ' - ' num2str(num)])
        %saveas(gcf, ['data/timelapse/figures/summedIm_' region '_' num2str(num) '.pdf'])
        
        scale = 250/1024; % um for 1 pixel
        if n == 3
            sc = 0.25;
        else sc = 0.5;
        end
        minSiz = sc*floor(1024/25)^2; % sc*10um in pixels (?? FELIX ??)
        
        % Detect on the concatenated video
        [cellMap1,nROI1] = cellDetectFun(cat(3,data1,data2),minSiz,0.5);

        % Plot all detected astros
        cmap = [0 0 0; hsv(max(nROI1,nROI1))];
        
        figure
        %imshow(gp_display_out_2p(1:nROI1,cellMap1,'off'))
        astrosGpeAll1 = zeros(size(gp_display_out_2p(1,cellMap1,'off')));

        for a = 1:nROI1
            astroGpe = gp_display_out_2p(a,cellMap1,'off');
            astrosGpeAll1 = astrosGpeAll1 + astroGpe*(a);   
        end

        figure
        subplot(2,2,1)
        colormap(cmap);
        image(astrosGpeAll1);
        title('All detected astrocytes')

        % Extract signals
        disp('Extracting signals')

        dt = 1/2;
        dSiz = 1;
        varMin = 0.1;
        [decharges1,dechNorm1,varList1] = extractIndivSignals(data1,cellMap1,dSiz,1:nROI1,varMin);
        [decharges2,dechNorm2,varList2] = extractIndivSignals(data2,cellMap1,dSiz,1:nROI1,varMin);
        totDur1 = length(decharges1(:,1));
        totDur2 = length(decharges2(:,1));
        
        % Define Global events
        dGlobEv1 = floor(totDur1/10);
        GlobalEventsTime1 = (1:dGlobEv1:totDur1)';
        evNumVec1 = 1:totDur1;
        for evN = 1:length(GlobalEventsTime1)
            evNumVec1((evN-1)*dGlobEv1+1:min(dGlobEv1*evN,totDur1))=evN;
        end
        
        dGlobEv2 = floor(totDur2/10);
        GlobalEventsTime2 = (1:dGlobEv2:totDur2)';
        evNumVec2 = 1:totDur2;
        for evN = 1:length(GlobalEventsTime2)
            evNumVec2((evN-1)*dGlobEv2+1:min(dGlobEv2*evN,totDur2))=evN;
        end
        
        % Plot the active astros
        nAstros1 = length(varList1);
        cmap = [0 0 0; hsv(nAstros1)];
        %imshow(gp_display_out_2p(varList1,cellMap1,'off'))
        astrosGpeAll1 = zeros(size(gp_display_out_2p(varList1(1),cellMap1,'off')));
        % Select the active astrocyte for event #ev
        for a = 1:nAstros1
            astroGpe = gp_display_out_2p(varList1(a),cellMap1,'off');
            astrosGpeAll1 = astrosGpeAll1 + astroGpe*(varList1(a));   
        end
        % Active astrocytes
        subplot(2,2,3)
        colormap(cmap);
        image(astrosGpeAll1);
        title('Active astros for t 1')
        
        subplot(2,2,4)
        %imshow(gp_display_out_2p(varList2,cellMap1,'off'))
        astrosGpeAll2 = zeros(size(gp_display_out_2p(varList2(1),cellMap1,'off')));
        % Select the active astrocyte for event #ev
        for a = 1:length(varList2)
            astroGpe = gp_display_out_2p(varList2(a),cellMap1,'off');
            astrosGpeAll2 = astrosGpeAll2 + astroGpe*(varList2(a));   
        end
        image(astrosGpeAll2);
        colormap(cmap);
        title('Active astros for t 2')
        
        subplot(2,2,2)
        image(astrosGpeAll1-astrosGpeAll2);
        colormap(cmap);
        title('Difference t1-t2')
        %saveas(gcf, ['data/timelapse/figures/detectedAstros_' region '_' num2str(num) '.pdf'])
        
        % Detect indiv events 
        disp('Detecting Individual events')

        dtev = 2;
        minDepth = 1.5;
        maxEvThresh = 1.5;
        Tev = 2.5;
        cfDown = 1.5;
        smoothParam = 7e-5;

        [meanFreq,realEvts1,diffToBaseAll,varL1] = ...
            segWithBaseCorr(decharges1(:,varList1),dt,dtev,minDepth,maxEvThresh,Tev,cfDown,...
            smoothParam,evNumVec1,GlobalEventsTime1);
        close all
        [meanFreq2,realEvts2,diffToBaseAll2,varL2] = ...
            segWithBaseCorr(decharges2(:,varList2),dt,dtev,minDepth,maxEvThresh,Tev,cfDown,...
            smoothParam,evNumVec2,GlobalEventsTime2);
        close all

        % Correlations
        disp('Doing correlations')
        dT = 0.01;
        [cr,corrCurve,corThresh] = pairCrossCorr(diffToBaseAll,varList1,dT);
        [cr,corrCurve2,corThresh] = pairCrossCorr(diffToBaseAll2,varList2,dT);

        figure
        plot(corThresh,corrCurve,'k','LineWidth',1.5)
        hold on 
        plot(corThresh,corrCurve2,'b','LineWidth',1.5)
        legend('t1','t2')
        xlabel('correlation threshold')
        ylabel('Proportion of astrocytes (%)')
        title('Proportion of astrocytes over correlation threshold')
        %saveas(gcf, ['data/timelapse/figures/corrCurves_' region '_' num2str(num) '.pdf'])

        % Compute AUC differencen and ks distance
        auc1 = trapz(corThresh, corrCurve);
        auc2 = trapz(corThresh, corrCurve2);
        diff_auc = abs(auc1-auc2);
        [h,p,ks2stat] = kstest2(corrCurve,corrCurve2);
        
        if strcmp(reg,' Cx Mtr ')
            aucCTX = [aucCTX,  diff_auc];
            ksCTX = [ksCTX,  ks2stat];
        elseif strcmp(reg,' CA1 ')
            aucCA1 = [aucCA1,  diff_auc];
            ksCA1 = [ksCA1,  ks2stat];
        end
        
        % Global events co activations and events paths
        disp('Doing paths')
        GloEvMat1 = globalEvtActiv(GlobalEventsTime1,realEvts1,varList1);
        GloEvMat2 = globalEvtActiv(GlobalEventsTime2,realEvts2,varList2);

        [paths,area,numActiv,evPathLength] =  activationPaths(GlobalEventsTime1,GloEvMat1,varList1,cellMap1);
        [paths2,area2,numActiv2,evPathLength2] =  activationPaths(GlobalEventsTime2,GloEvMat2,varList2,cellMap1);
        close all
        if strcmp(reg,' Cx Mtr ')
            areaCTX = [areaCTX,  [mean(area); mean(area2)]];
            numActivCTX = [numActivCTX,  [mean(numActiv); mean(numActiv2)]];
            evPathLengthCTX = [evPathLengthCTX, [mean(evPathLength); mean(evPathLength2)]];
        elseif strcmp(reg,' CA1 ')
            areaCA1 = [areaCA1,  [mean(area); mean(area2)]];
            numActivCA1 = [numActivCA1,  [mean(numActiv); mean(numActiv2)]];
            evPathLengthCA1 = [evPathLengthCA1, [mean(evPathLength); mean(evPathLength2)]];
        end
        
        % Graphs
        disp('Doing graphs')
        T_HC = 10;
        orientation = 'no';
        [netGraph,meanDeg,varDeg,numHighConnNodes,propHighConnNodes,conMat] = ...
            activGraph(GlobalEventsTime1,GloEvMat1,cellMap1,varList1,T_HC,orientation,10/scale);
        saveas(gcf, ['data/timelapse/figures/graph_' region '_' num2str(num) '_t1.pdf'])
        
        [netGraph2,meanDeg2,varDeg2,numHighConnNodes2,propHighConnNodes2,conMat2] = ...
            activGraph(GlobalEventsTime2,GloEvMat2,cellMap1,varList2,T_HC,orientation,10/scale);
        %saveas(gcf, ['data/timelapse/figures/graph_' region '_' num2str(num) '_t2.pdf'])
        if strcmp(reg,' Cx Mtr ')
            propHCCTX = [propHCCTX, [propHighConnNodes; propHighConnNodes2]];
            numHCCTX = [numHCCTX, [numHighConnNodes; numHighConnNodes2]];
            meanDegCTX = [meanDegCTX, [meanDeg; meanDeg2]];
            varDegCTX = [varDegCTX, [varDeg; varDeg2]];
        elseif strcmp(reg,' CA1 ')
            propHCCA1 = [propHCCA1, [propHighConnNodes; propHighConnNodes2]];
            numHCCA1 = [numHCCA1, [numHighConnNodes; numHighConnNodes2]];
            meanDegCA1 = [meanDegCA1, [meanDeg; meanDeg2]];
            varDegCA1 = [varDegCA1, [varDeg; varDeg2]];
        end
        
        % Connectivity matrix
        disp('Doing conn mats')
        clim=[0,20];
        figure
        imagesc(conMat,clim)
        title('t1')
        colormap hot
        colorbar        
        %saveas(gcf, ['data/timelapse/figures/conMats_' region '_' num2str(num) '_t1.pdf'])
        
        figure
        imagesc(conMat2,clim)
        title('t2')
        colormap hot
        colorbar
        %saveas(gcf, ['data/timelapse/figures/conMats_' region '_' num2str(num) '_t2.pdf'])
        ssimMats = ssim(conMat,conMat2);
        if strcmp(reg,' Cx Mtr ')
            ssimCTX = [ssimCTX,  ssimMats];
        elseif strcmp(reg,' CA1 ')
            ssimCA1 = [ssimCA1,  ssimMats];
        end
    end 
end
close all
%% Plots the results t1 vs t2
cmap = hsv;

CTXcol = 10;
HPCcol = 50;

% Graph measures t1 vs t2
% Prop HC
OvecCA1 = ones(length(propHCCA1),1);
figure
subplot(3,2,1)
plot(OvecCA1,propHCCA1(1,:),'.','MarkerSize',20,'color',cmap(HPCcol,:))
hold on 
plot(2*OvecCA1,propHCCA1(2,:),'.','MarkerSize',20,'color',cmap(HPCcol,:))
plot([OvecCA1 2*OvecCA1]', [propHCCA1(1,:); propHCCA1(2,:)],'-','color',cmap(HPCcol,:))

OvecCTX = ones(length(propHCCTX),1);
plot(OvecCTX,propHCCTX(1,:),'.','MarkerSize',20,'color',cmap(CTXcol,:))
hold on 
plot(2*OvecCTX,propHCCTX(2,:),'.','MarkerSize',20,'color',cmap(CTXcol,:))
plot([OvecCTX 2*OvecCTX]', [propHCCTX(1,:); propHCCTX(2,:)],'-','color',cmap(CTXcol,:))
title('Proportion of HC nodes')
xlim([0,3])

% Num HC
subplot(3,2,2)
plot(OvecCA1,numHCCA1(1,:),'.','MarkerSize',20,'color',cmap(HPCcol,:))
hold on 
plot(2*OvecCA1,numHCCA1(2,:),'.','MarkerSize',20,'color',cmap(HPCcol,:))
plot([OvecCA1 2*OvecCA1]', [numHCCA1(1,:); numHCCA1(2,:)],'-','color',cmap(HPCcol,:))

plot(OvecCTX,numHCCTX(1,:),'.','MarkerSize',20,'color',cmap(CTXcol,:))
hold on 
plot(2*OvecCTX,numHCCTX(2,:),'.','MarkerSize',20,'color',cmap(CTXcol,:))
plot([OvecCTX 2*OvecCTX]', [numHCCTX(1,:); numHCCTX(2,:)],'-','color',cmap(CTXcol,:))
title('Number of HC nodes')
xlim([0,3])

% Mean node degree
subplot(3,2,3)
plot(OvecCA1,meanDegCA1(1,:),'.','MarkerSize',20,'color',cmap(HPCcol,:))
hold on 
plot(2*OvecCA1,meanDegCA1(2,:),'.','MarkerSize',20,'color',cmap(HPCcol,:))
plot([OvecCA1 2*OvecCA1]', [meanDegCA1(1,:); meanDegCA1(2,:)],'-','color',cmap(HPCcol,:))

plot(OvecCTX,meanDegCTX(1,:),'.','MarkerSize',20,'color',cmap(CTXcol,:))
hold on 
plot(2*OvecCTX,meanDegCTX(2,:),'.','MarkerSize',20,'color',cmap(CTXcol,:))
plot([OvecCTX 2*OvecCTX]', [meanDegCTX(1,:); meanDegCTX(2,:)],'-','color',cmap(CTXcol,:))
title('Mean node degree')
xlim([0,3])

% Path info
% Evt area
subplot(3,2,4)
plot(OvecCA1,areaCA1(1,:),'.','MarkerSize',20,'color',cmap(HPCcol,:))
hold on 
plot(2*OvecCA1,areaCA1(2,:),'.','MarkerSize',20,'color',cmap(HPCcol,:))
plot([OvecCA1 2*OvecCA1]', [areaCA1(1,:); areaCA1(2,:)],'-','color',cmap(HPCcol,:))

plot(OvecCTX,areaCTX(1,:),'.','MarkerSize',20,'color',cmap(CTXcol,:))
hold on 
plot(2*OvecCTX,areaCTX(2,:),'.','MarkerSize',20,'color',cmap(CTXcol,:))
plot([OvecCTX 2*OvecCTX]', [areaCTX(1,:); areaCTX(2,:)],'-','color',cmap(CTXcol,:))
title('Mean active area')
xlim([0,3])

% Path length
subplot(3,2,5)
plot(OvecCA1,evPathLengthCA1(1,:),'.','MarkerSize',20,'color',cmap(HPCcol,:))
hold on 
plot(2*OvecCA1,evPathLengthCA1(2,:),'.','MarkerSize',20,'color',cmap(HPCcol,:))
plot([OvecCA1 2*OvecCA1]', [evPathLengthCA1(1,:); evPathLengthCA1(2,:)],'-','color',cmap(HPCcol,:))

plot(OvecCTX,evPathLengthCTX(1,:),'.','MarkerSize',20,'color',cmap(CTXcol,:))
hold on 
plot(2*OvecCTX,evPathLengthCTX(2,:),'.','MarkerSize',20,'color',cmap(CTXcol,:))
plot([OvecCTX 2*OvecCTX]', [evPathLengthCTX(1,:); evPathLengthCTX(2,:)],'-','color',cmap(CTXcol,:))
title('Mean path length')
xlim([0,3])

% Num activ
subplot(3,2,6)
plot(OvecCA1,numActivCA1(1,:),'.','MarkerSize',20,'color',cmap(HPCcol,:))
hold on 
plot(2*OvecCA1,numActivCA1(2,:),'.','MarkerSize',20,'color',cmap(HPCcol,:))
plot([OvecCA1 2*OvecCA1]', [numActivCA1(1,:); numActivCA1(2,:)],'-','color',cmap(HPCcol,:))

plot(OvecCTX,numActivCTX(1,:),'.','MarkerSize',20,'color',cmap(CTXcol,:))
hold on 
plot(2*OvecCTX,numActivCTX(2,:),'.','MarkerSize',20,'color',cmap(CTXcol,:))
plot([OvecCTX 2*OvecCTX]', [numActivCTX(1,:); numActivCTX(2,:)],'-','color',cmap(CTXcol,:))
title('Number active astrocytes')
xlim([0,3])

%% Diff measures
% ssim conMats
figure
subplot(1,2,1)
plot(OvecCA1,ssimCA1,'*','MarkerSize',10,'color',cmap(HPCcol,:))
hold on
plot(OvecCTX,ssimCTX,'*','MarkerSize',10,'color',cmap(CTXcol,:))
boxplot([ssimCA1, ssimCTX])
title('SSIM between connectivity matrices')
ylim([0,1])
% AUC diff corrCurves
subplot(1,2,2)
plot(OvecCA1,aucCA1,'*','MarkerSize',10,'color',cmap(HPCcol,:))
hold on
plot(OvecCTX,aucCTX,'*','MarkerSize',10,'color',cmap(CTXcol,:))
boxplot([aucCA1, aucCTX])
title('AUC difference')
ylim([0,1])