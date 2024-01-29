% This script processes shows an example of how to process several
% recorded sessions (in this example the video data is already converted to
% .mat format but the script works for any matlab supported video format, 
% you should just adapt step 0 below for the loading part).
% files and concatenates the results (steps 1 and 2) into structure arrays 
% that can be used for the different types of data analysis proposed (steps 3A-B-C)
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% L. Zonca, Jan. 2024
% ----------------------------------------------------------------------- %
% Choose Data path
path = 'C:\Users\louzo\data_path_example\';
% Load sessions names
load([path 'allSessionsNames.mat']);

% Choose to make plots or not
printPlots = 'no'; % 'yes'; % 

% Initialize structure
numOfSessions = 5;
allDATA = cell(numOfSessions,1); % This cell array will contain all the results

% ---- Fixed Parameters --- %
cmap = hsv;
% Cell detection
scale = 250/1024; % in um for 1 pixel
minSiz = floor(1024/25)^2; % e.g. 10um in pixels
% Signal extraction
dt = 1/2; % in seconds
dSiz = 1;
varMin = 0.1;
% Indiv. event detection
dtev = 2;
minDepth = 3;
maxEvThresh = 3;
Tev = 2.5;
cfDown = 1.5;
smoothParam = 7e-5;

% Highly connected threshold
T_HC = 12; % 60% de connectivité 
orientation = 'no';

%% Run the processing
for sessionIdx = 1:numOfSessions
    disp(['Processing session #' num2str(sessionIdx)])
    session = allSessionNames{sessionIdx};

    % 0 - Load data
    dataLoad = load([path '\' session '.mat']); % Replace with data location
    data = dataCorr.data;
    % 1 - Detect cells
    [cellMap,nROI] = cellDetectFun(data,minSiz,0.5);
    close all
    % 2 - Extract signals
    [decharges,dechNorm,varList] = extractIndivSignals(data,cellMap,dSiz,1:nROI,varMin);
    totDur = length(decharges(:,1));
    timeVec = 1:dt:totDur*dt;
    nAstros = length(varList);
    close all

    % 2.1 - Define Global events
    dGlobEv = floor(totDur/10);
    GlobalEventsTime = (1:dGlobEv:totDur)';
    evNumVec = 1:totDur;
    for evN = 1:length(GlobalEventsTime)
        evNumVec((evN-1)*dGlobEv+1:min(dGlobEv*evN,totDur))=evN;
    end

    % 3 - Detect indiv events 
    [meanFreq,realEvts,diffToBaseAll,varList2] = ...
        segWithBaseCorr(decharges(:,varList),...
        dt,dtev,minDepth,maxEvThresh,Tev,cfDown,...
        smoothParam,evNumVec,GlobalEventsTime);
    close all

    % 4 - Make correlations curves
    dT = 0.01;
    [cr,corrCurve,corThresh] = pairCrossCorr(diffToBaseAll,varList,dT);

    % 5 - Global events co-activations and events paths
    GloEvMat = globalEvtActiv(GlobalEventsTime,realEvts,varList);
    [paths,area,numActiv,evPathLength,pathStepsAll,pathStepsTimes] =  ...
         activationPaths(GlobalEventsTime,GloEvMat,varList,cellMap);

    % 6 - Connectivity Graphs
    [netGraph,meanDeg,varDeg,numHighConnNodes,propHighConnNodes, conMat] = ...
        activGraph(GlobalEventsTime,GloEvMat,cellMap,varList,T_HC,orientation,10/scale);
    close all
    % Add results in global structure array
    sessionResults = struct('name',session{sess},...
        'timeVec',timeVec,...
        'GlobalEventsTime',GlobalEventsTime,...
        'allActiveAstros',varList,...
        'nAstros',nAstros,...
        'meanEvtFreq',meanFreq,...
        'evtBegin',{realEvts(:,1)},...
        'evtPeak',{realEvts(:,2)},...
        'evtEnd',{realEvts(:,3)},...
        'evtGlobal',{realEvts(:,4)},...
        'evtAmpli',{realEvts(:,5)},...
        'evtSubPeaks',{realEvts(:,6)},...
        'evtSubPeaksFreq',{realEvts(:,7)},...
        'corrCurve',corrCurve,...
        'activeAreas',area,...
        'activeAstrosEvt',numActiv,...
        'pathLengths',evPathLength,...
        'meanNodeDeg',meanDeg,...
        'varNodeDeg',varDeg,...
        'numHighConnNodes',numHighConnNodes,...
        'propHighConnNodes',propHighConnNodes,...
        'conMat',conMat,...
        'pathStepsTimes',{pathStepsTimes},...
        'pathSteps',{pathStepsAll});
    allDATA{sessionIdx,sess} = sessionResults;
    disp([session{sess} ' loaded (session #' num2str(sess) '/' num2str(length(session)) ')'])

end


save('allDATA.mat','allDATA')