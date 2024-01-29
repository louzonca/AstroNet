function [netGraph,meanDeg,varDeg,numHighConnNodes,propHighConnNodes,conMat] = ...
    activGraph(GlobalEventsTime,GloEvMat,cellMap,varList,T_HC,orientation,scale)
% This function builds the activation graph for one recording session
% ----------------------------------------------------------------------- %
% *** Inputs ***
% * GlobalEventsTime * is a vector containing the time points (in frames)
% delimiting the global events
% * GloEvMat * is a cell array containing the list of active ROIs and their 
% peak time for each global event.
% * cellMap * is the map of detected ROIs
% * varList * is the list of active ROIs
% * T_HC * is the threshold value for highly connected nodes
% * orientation * is a string = 'yes' for an oriented graph and 'no' for a
% non oriented graph
% * scale * is the spatial scale in pixel
% ----------------------------------------------------------------------- %
% *** Outputs ***
% * netGraph * is the graph
% * meanDeg * is the mean node degree in *netGraph*
% * varDeg * is the standard deviation of node degrees in *netGraph*
% * numHighConnNodes * is the number of highly connected nodes i.e. with at
% least one edge of weight >= *T_HC* 
% * propHighConnNodes * is the proportion (in %) of highly connected nodes
% * conMat * is the connectivity matrix between the ROIs (directed if 
% *orientation* = 'yes')
% ----------------------------------------------------------------------- %
% L. Zonca, Jan. 2022
% ----------------------------------------------------------------------- %
nAstros = length(varList);
% Centroids positions
s = regionprops(cellMap, {'Centroid', 'PixelIdxList'});
centroid = zeros(max(varList),2);
for k = 1:max(varList)
   centroid(k,:) = s(k).Centroid;
end
centroid(setdiff(1:max(varList),varList),:) = 0;

% --- Option 1: Non oriented graph ---
if strcmp(orientation,'no')
    % --- Build non oriented connectivity matrix ---
    conMat = zeros(varList(end));
    for ev = 1:length(GlobalEventsTime)
        try
            conPairs = nchoosek(GloEvMat{ev,1},2);
            conPairsMat = zeros(varList(end));
            conPairsMat(conPairs(:,1),conPairs(:,2)) = 1;
            conMat = conMat + conPairsMat + conPairsMat';
        catch
            conPairsMat = zeros(varList(end));
            conPairsMat(GloEvMat{ev,1},GloEvMat{ev,1}) = 0;
            conMat = conMat + conPairsMat + conPairsMat';
        end
    end
    %  --- Make the non oriented graph ---
    G =  graph(conMat,'omitselfloops');
    netGraph = G;
    %%%%%%%%%%%%%% Display %%%%%%%%%%%%%%%
    % if strcmp(printPlots,'yes')
        figure 
        % ax1 = subplot(1,2,1);
        % Colormap
        cmapAs = [0 0 0; hsv(nAstros)];
        cmap = hot(2*(length(GlobalEventsTime)));%hot(max(G.Edges.Weight));
        colTot = [cmapAs; cmap];

        % Plot the active ROIs --- Non oriented graph ---
        astrosGpeAll = zeros(size(gp_display_out_2p(varList(1),cellMap,'off')));
        for a = 1:length(varList)
            astroGpe = gp_display_out_2p(varList(a),cellMap,'off');
            astrosGpeAll = astrosGpeAll + astroGpe*varList(a);
        end
        colormap(colTot);
        image(astrosGpeAll);
        hold on

        % Plot the graph --- Non oriented graph ---
        LWidth = 10*G.Edges.Weight/max(G.Edges.Weight);
        EdgCol = colTot(nAstros+1+G.Edges.Weight,:);
        plot(G,'LineWidth',LWidth,'EdgeColor',EdgCol,'Xdata',centroid(1:varList(end),1),'YData',centroid(1:varList(end),2));
    % end
    %%%%%%%%%%%%% End of display %%%%%%%%%%%%%%

    % Degree of each graph node --- Non oriented graph ---
    D = degree(G);
    meanDeg = mean(D(varList));
    varDeg = std(D(varList));
%     if strcmp(printPlots,'yes')
        colormap(gca,colTot);
        cb = colorbar('location','southoutside','Ticks',nAstros+2:length(colTot),...
                 'TickLabels',1:length(colTot)-nAstros-1);
        set(cb, 'ylim', [nAstros+2 length(colTot)])
%     end
    % Find nodes with at least one edge with weight > T_HC --- Non oriented graph ---
    highConnNodes = [];
    numHighConnNodes = 0;
    for n = 1:length(conMat(1,:))
        if max(conMat(n,:)) >= T_HC
            highConnNodes = [highConnNodes; n];
            numHighConnNodes = numHighConnNodes +1;
        end
    end
    propHighConnNodes = numHighConnNodes/numel(varList)*100;
    % Indicate proportion of high connected nodes on graph
    text(15,15,...
        ['Highly connected nodes: ' num2str(numHighConnNodes/numel(varList)*100,'%4.2f') '%'],...
        'Color','white','FontSize',12)
    % Add spatial scale
    plot([size(astrosGpeAll,1)-50 size(astrosGpeAll,1)-50+scale],...
        [size(astrosGpeAll,2)-50 size(astrosGpeAll,2)-50],'color','w','LineWidth',4)
    text(size(astrosGpeAll,1)-50-scale,size(astrosGpeAll,2)-100,'10 um','Color','white','FontSize',12)

% --- Option 2: Oriented graph ---
elseif strcmp(orientation,'yes')
    % --- Build directed connectivity matrix for oriented graph ---
    dirConMat = zeros(varList(end));
    for ev = 1:length(GlobalEventsTime)
        try
            conPairs = nchoosek(GloEvMat{ev,1},2);       
            for pn = 1:length(conPairs(:,1))
                conPairsMat = zeros(varList(end));
                if GloEvMat{ev,2}(find(GloEvMat{ev,1} == conPairs(pn,1)))>=...
                        GloEvMat{ev,2}(find(GloEvMat{ev,1} == conPairs(pn,2)))
                        conPairsMat(conPairs(pn,1),conPairs(pn,2)) = 1;
                else conPairsMat(conPairs(pn,2),conPairs(pn,1)) = 1;
                end
                dirConMat = dirConMat + conPairsMat;
            end
        catch
            conPairsMat = zeros(varList(end));
            conPairsMat(GloEvMat{ev,1},GloEvMat{ev,1}) = 0;
            dirConMat = dirConMat + conPairsMat;
        end
    end

    % --- Make the oriented graph ---
    Gdir =  digraph(dirConMat,'omitselfloops');
    netGraph = Gdir;
    %%%%%%%%%%%%%% Display oriented graph %%%%%%%%%%%%%%%
    % if strcmp(printPlots,'yes')
        % ax2 = subplot(1,2,2);
        % Colormap --- Oriented graph ---
        cmapAs = [0 0 0; hsv(nAstros)];
        cmap = hot(max(Gdir.Edges.Weight));
        colTotDir = [cmapAs; cmap];

        % Plot the active ROIs --- Oriented graph ---
        astrosGpeAll = zeros(size(gp_display_out_2p(varList(1),cellMap,'off')));
        for a = 1:length(varList)
            astroGpe = gp_display_out_2p(varList(a),cellMap,'off');
            astrosGpeAll = astrosGpeAll + astroGpe*varList(a);
        end
        colormap(colTotDir);
        image(astrosGpeAll);
        hold on

        LWidth = 5*Gdir.Edges.Weight/max(Gdir.Edges.Weight);
        EdgCol = colTotDir(nAstros+1+Gdir.Edges.Weight,:);
        p = plot(Gdir,'LineWidth',LWidth,'EdgeColor',EdgCol,'Xdata',centroid(1:varList(end),1),'YData',centroid(1:varList(end),2));
        p.ArrowSize = 15;
        colormap(gca,colTotDir);
        cbDir = colorbar('location','southoutside','Ticks',nAstros+2:length(colTotDir),...
                 'TickLabels',1:length(colTotDir)-nAstros-1);
        set(cbDir, 'ylim', [nAstros+2 length(colTotDir)])
    % end
    %%%%%%%%%%% End of display oriented graph %%%%%%%%%%
    
    % Degree of each graph node --- Ooriented graph ---
    Din = indegree(Gdir);
    Dout = outdegree(Gdir);
    meanDegIn = mean(Din(varList));
    varDegIn = std(Din(varList));
    meanDegOut = mean(Dout(varList));
    varDegOut = std(Dout(varList));
    meanDeg = [meanDegIn, meanDegOut];
    varDeg = [varDegIn, varDegOut];
    % Find nodes with at least one edge with weight > T_HC --- Non oriented graph ---
    T_HC_dir = ceil(T_HC/2);
    highConnNodes = [];
    numHighConnNodes = 0;
    for n = 1:length(dirConMat(1,:))
        if max(abs(dirConMat(n,:))) >= T_HC_dir
            highConnNodes = [highConnNodes; n];
            numHighConnNodes = numHighConnNodes +1;
        end
    end
    propHighConnNodes = numHighConnNodes/numel(varList)*100;
    % Indicate proportion of high connected nodes on graph
    text(15,15,...
        ['Highly connected nodes: ' num2str(numHighConnNodes/numel(varList)*100,'%4.2f') '%'],...
        'Color','white','FontSize',12)
    % Add spatial scale
    plot([size(astrosGpeAll,1)-50 size(astrosGpeAll,1)-50+scale],...
        [size(astrosGpeAll,2)-50 size(astrosGpeAll,2)-50],'color','w','LineWidth',4)
    text(size(astrosGpeAll,1)-50-scale,size(astrosGpeAll,2)-100,'10 um','Color','white','FontSize',12)
    conMat = dirConMat;
    
else warning('orientation must be yes or no')
end