function [paths,area,numActiv,evPathLength,pathStepsAll,pathStepsTimes] = ...
    activationPaths(GlobalEventsTime,GloEvMat,varList,cellMap)
% This function builds the activation paths for all global events in one
% recording session
% ----------------------------------------------------------------------- %
% *** Inputs ***
% * GlobalEventsTime * is a vector containing the time points (in frames)
% delimiting the global events
% * GloEvMat * is a cell array containing the list of active ROIs and their 
% peak time for each global event.
% * varList * is the list of active ROIs
% * cellMap * is an image giving the position and label# of all ROIs
% ----------------------------------------------------------------------- %
% *** Outputs ***
% * paths * is a cell array containing the ordered list of ROI centroid
% positions along the activation path
% * area * is a vector containing the area (in square pixels) of the convex 
% hull of the active ROIs in each global event
% * numActiv * is a vector containing the number of active ROIs per global event 
% * evPathLength * is a vector containing the total length (in pixels) of 
% the activation path for each global event.
% * pathStepsAll * is a cell array containing a vector of consecutive path
% steps length for each path
% * pathStepsTimes * is a cell array containing a vector of consecutive path
% steps durations for each path
% ----------------------------------------------------------------------- %
% L. Zonca, Jan. 2022
% ----------------------------------------------------------------------- %
nAstros = length(varList);
cmap = [0 0 0; hsv(nAstros)];

area = zeros(size(GlobalEventsTime));
numActiv = zeros(size(GlobalEventsTime));
evPathLength = zeros(size(GlobalEventsTime));
paths = cell(length(GlobalEventsTime),1);
pathStepsAll = cell(length(GlobalEventsTime),1);
pathStepsTimes = cell(length(GlobalEventsTime),1);

% Centroids positions
s = regionprops(cellMap, {'Centroid', 'PixelIdxList'});
centroid = zeros(max(varList),2);
for k = 1:max(varList)
   centroid(k,:) = s(k).Centroid;
end
centroid(setdiff(1:max(varList),varList),:) = 0;

for ev = 1:length(GlobalEventsTime)
    [~,order] = sort(GloEvMat{ev,2});
    astrosGpeAll = zeros(size(gp_display_out_2p(varList(1),cellMap,'off')));
    % Select the active astrocyte for event #ev
    for a = 1:length(GloEvMat{ev,1})
        astroGpe = gp_display_out_2p(GloEvMat{ev,1}(a),cellMap,'off');
        astrosGpeAll = astrosGpeAll + astroGpe*(GloEvMat{ev,1}(a));   
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(sum(astrosGpeAll))>0 && length(GloEvMat{ev,1}) >1
        figure
        % Active astrocytes
        colormap(cmap);
        image(astrosGpeAll);
        % Activation paths
        pathLen = 0;
        hold on
        path = zeros(length(GloEvMat{ev,1}),2);
        pathSteps = zeros(length(GloEvMat{ev,1})-1,1);
        pathStepsT = zeros(length(GloEvMat{ev,1})-1,1); 
        for a = 1:length(GloEvMat{ev,1})-1
            asNum = order(a);%find(order==a,1);
            asNum1 = order(a+1);%find(order==a+1,1);
            x = centroid(GloEvMat{ev,1}(asNum),:);
            path(a,:) = x;
            y = centroid(GloEvMat{ev,1}(asNum1),:); 
            dp = (y-x)/2;
            cmap2 = [0 0 0; hot(length(GloEvMat{ev,2})+1)];
            plot([x(1) y(1)],[x(2) y(2)],'LineWidth',2.,...
                'color',cmap2(a+1,:))
            % Plot the arrows
            q = quiver(x(1),x(2),dp(1),dp(2),0,'color',cmap2(a+1,:),'LineWidth',2.,'MarkerSize',20);          
            % Path length
            pathSteps(a) = sqrt((centroid(GloEvMat{ev,1}(asNum),1)-centroid(GloEvMat{ev,1}(asNum1),1))^2 +...
                (centroid(GloEvMat{ev,1}(asNum),2)-centroid(GloEvMat{ev,1}(asNum1),2))^2);
            pathLen = pathLen + pathSteps(a);
            % Path steps times
            pathStepsT(a) = (GloEvMat{ev,2}(asNum1)-GloEvMat{ev,2}(asNum));
            %get the data from regular quiver
            U = q.UData;
            V = q.VData;
            X = q.XData;
            Y = q.YData;
            % Overlap arrows with annotation
            headLength = 12.5;
            headWidth = 8;
            for ii = 1:length(X)
                for ij = 1:length(X)
                    set(gca,'YDir','normal')
                    ah = annotation('arrow','color',cmap2(a+1,:),...
                        'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
                    set(ah,'parent',gca);
                    set(ah,'position',[X(ii,ij) Y(ii,ij) U(ii,ij) V(ii,ij)]);
                end
            end
        end
        path(a+1,:) = y;
        % Add circle around the first and last astrocyte of the path
        % first
        asNumFirst = order(1);
        pos = [centroid(GloEvMat{ev,1}(asNumFirst),1)-11 ...
            centroid(GloEvMat{ev,1}(asNumFirst),2)-11 22 22];
        rectangle('Position',pos,'Curvature',[1 1],'EdgeColor','g','LineWidth',2)
        % Last
        asNumLast = order(length(GloEvMat{ev,1}));
        pos = [centroid(GloEvMat{ev,1}(asNumLast),1)-11 ...
            centroid(GloEvMat{ev,1}(asNumLast),2)-11 22 22];
        rectangle('Position',pos,'Curvature',[1 1],'EdgeColor','r','LineWidth',2)
           
        evPathLength(ev) = pathLen;
        % Add convex hull of active astrocytes
        try
            [borders,area(ev)] = convhull(centroid(GloEvMat{ev,1},1),centroid(GloEvMat{ev,1},2));
                plot(centroid(GloEvMat{ev,1}(borders),1),centroid(GloEvMat{ev,1}(borders),2),...
                    'w--','LineWidth',1.5)
        catch
            area(ev) = 0;
        end
        % Number of active astrocyte in event
        numActiv(ev) = length(GloEvMat{ev,1});  
        title(['Event ' num2str(ev) ' active area = ' num2str(area(ev)) ' path length ' num2str(evPathLength(ev))]) 
        paths{ev} = path;
        pathStepsAll{ev} = pathSteps;
        pathStepsTimes{ev} = pathStepsT;
    end   
end