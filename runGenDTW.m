function runGenDTW(measure)

close all
addPath

if nargin == 0
    measure = 'overlap';
end

%% Global variables
global DATATYPE;        % dataset: 'chalearn' 'random'
global VISUALIZE;       % flag indicating whether or not to visualize the data groups using kmeans clustering
global FRAMECOMPRESS;   % flag indicating whether or not to perform frame compression 
global ALLLEARNING;     % flag indicating whether or not to perform frame compression 
global KMEANSDTWv;      % possible versions of the k-means DTW algorithm
global DISTANCES;       % possible distances: 'dtw' 'euclidean'
global SAMPLING;        % sampling type
global S;               % Scores from the GA training
global CACHE;           % Cache for the GA
% global GENRESULTS;      % Population at each generation and results
global PREDICTIONS;     % Predictions reached: [startFrame,endFrame,predLabel]
global BASELINE;        % Baseline to be used in the GA
global COORDS;          % type of coordinates
global STATE;           % save state of previous GA run
global OPTIONS;         % Options of the GA run
global PERCENTDATA;     % percent of learning data to consider
global NORMTYPE;        % Normalization type 'neck' 'xyzangles'
global MEDIANTYPE;      % Type of median models 
global JOINTS;          % Selected joints
global NAT;             % Natalia's descriptor

DATATYPE = 'chalearn2014';
NORMTYPE = 'neck';
COORDS = 'pixel';
VISUALIZE = false;
FRAMECOMPRESS = true;
ALLLEARNING = true;
KMEANSDTWv = {'v1_0','v1_1','v1_2','v1_3','v2_0_0','v2_0_1','v2_1_0',...
    'v2_1_1','v2_2_0','v2_2_1','v2_3_0','v2_3_1'};
DISTANCES = {'dtwCost','euclidean'};
SAMPLING = {'random','labels','segments'};
S = [];
BASELINE = {'random','deriv','fixed'};
STATE = [];
OPTIONS = [];
PERCENTDATA = 100;
MEDIANTYPE = {'direct','KNN','DCSR'};
JOINTS = [4 6 7 8 10 11 12];%1:20;%  %[8,12] hands
NAT = 0;

%% parameters data structures
nrsamples = 100;        % number of random samples
nseqs = 0;              % number of batches for training and validation, respectively. '0': All ; [nTrain nVal]
nframesSeg = 0;         % Fixed number of frames for subgesturing. '0' means no fixed subgesturing
sampling = SAMPLING{2}; % sampling type: 'random' 'label' 'segments'
noise = false;          % flag indicating whether or not consider noise for the test sequence
secsBatch = 60;         % reference seconds for the test sequence
nSampGest = 0;          % Number of samples per gesture for the test sequence   
mType = MEDIANTYPE{1};  % Type of median models to consider

% genetic temporal clustering parameters
params.versions = ...
    {KMEANSDTWv{5}};        % versions of the k-means DTW algorithm to execute';
params.dist = DISTANCES{1}; % distance metric
params.k0 = 3;              % initial data clusters for subgesturing
params.nmin = 5;            % minimum subsequence width
params.nmax = 25;           % maximum subsequence width
params.N = 500;             % Number of segments to split the learning sequence
params.nThreshs = 100;      % Number of thresholds for testing
params.D = [];              % Dissimilarity matrix
params.bestThs = [];        % Thresholds learnt on training
params.vectorized = 'on';   % vectorize the GA
params.population = 10;      % population of the GA
params.generations = 1000;  % number of generations of the GA
params.scoreMeasure ...
    = measure;            % Score Measure: 'overlap' or 'levenshtein'
params.Baseline = ...
    BASELINE{2};            % Baseline for the GA
params.threshMov = 3;       % maximum number of low movement frames
params.thMinMov = 1.3;      % Minimum portion of movement
params.drawMovSkels = false;% flag indicating whether to draw skels
params.thMutCross = 0.6;    % threshold of GA mutation and crossover
params.scale = 0.5;         % scale parameter for Gaussian mutation
params.shrink = 0.75;       % shrink parameter for Gaussian mutation
params.probSeg = 0.2;       % probability of eliminate/change a segment
params.maxWlen = 1000;      % maximum DTW cost matrix length to detect the start-end
if strcmp(mType,'KNN')
	params.k = 3;
else
	params.k = 0;           % current k to evaluate for the K-Nearest Neighbour models
end

CACHE.ind = int32(zeros(... % Cache with the populations
    params.population*100,params.N*2+1,'int32'));
if strcmp(measure,'overlap')
    CACHE.eval = zeros(...   % Cache with the scores of each individual
        params.population*100,1);
elseif strcmp(measure,'levenshtein')
    CACHE.eval = inf*ones(...% Cache with the scores of each individual
        params.population*100,1);
end
CACHE.pos = int32(1);       % Index positions
% GENRESULTS.P = cell(1,params.generations);
% GENRESULTS.eval = cell(1,params.generations);
PREDICTIONS = [];

%% Prepare training data depending on the chosen option and parameters
% Load data:
%     if nframesSeg is 0, then initial segmentation is generated from the skeleton labels
[X,Y,Xtest,Ytest] = prepareData(nrsamples,nseqs,nframesSeg,params.k0);
% display('Press a key to continue...');
% pause();

%% Compute initial segmentation from motion
[seg0,fr_fixed,params] = computeMagnitudes(X{1},params);
% [~,~,Xtrain,I,~,segTrain,~] = getDataSegments(X,Y,params.N,params.k0,params.nmin,params.nmax);

%% Prepare training data depending on the chosen option and parameters
DATATYPE = 'chalearn2014';
NORMTYPE = 'none';
COORDS = 'world';
NAT = 3;
% Load data:
%     if nframesSeg is 0, then initial segmentation is generated from the skeleton labels
[X,Y,Xtest,Ytest] = prepareData(nrsamples,nseqs,nframesSeg,params.k0);
% X = setDerivatives(X);
% display('Press a key to continue...');
% pause();

%% Obtain all samples grouped by gestures
Xtrain_l = getGroupedGestures(X,Y,1);
% Xval_l = getGroupedGestures(X,Y,2);
%Xtrain_l = getGroupedGestures(X,Y,0);

%% Get median models from training/learning data
fprintf('Computing mean gesture Models from training for the m=%d distinct gestures ...\n',length(Xtrain_l)-1);
% profile -memory on
params.M = getMedianModels(Xtrain_l,length(Xtrain_l)-1,mType,false);
% profreport
display('Done!');

%% Generate development sequences
% l = [24 78 150];    % 78 (more samples for each gesture when k=3);
l = [];
[Xdev,Ydev] = getDevSequences(X,Y,l,noise,secsBatch,nSampGest);
% Xdev{1} = Xdev{1}(1:round(length(Xdev{1})/4),:);
% Ydev{1}.Lfr = Ydev{1}.Lfr(1:round(length(Ydev{1}.Lfr)/4));
% Ydev{1}.L = unique(Ydev{1}.Lfr);
% Ydev{2}.Lfr = [Ydev{2}.Lfr(1002:1059) Ydev{2}.Lfr(280:324) ...
%     Ydev{2}.Lfr(233:279) Ydev{2}.Lfr(756:835) Ydev{2}.Lfr(964:1001)];
% Xdev{2} = [Xdev{2}(1002:1059,:); Xdev{2}(280:324,:); ...
%     Xdev{2}(233:279,:); Xdev{2}(756:835,:); Xdev{2}(964:1001,:)];
% Ydev{2}.L = unique(Ydev{2}.Lfr);

% Xdev = cell(1,2); Ydev = cell(1,2); 
% Ytrain_l = cell(1,length(Xtrain_l)-1); Yval_l = cell(1,length(Xval_l)-1);
% for i = 1:length(Xtrain_l)-1
%     Xtrain_l{i} = toMat(Xtrain_l{i});
%     Ytrain_l{i}.Lfr = i*ones(1,length(Xtrain_l{i}));
%     Xval_l{i} = toMat(Xval_l{i});
%     Yval_l{i}.Lfr = i*ones(1,length(Xval_l{i}));
% end
% Xtrain_l(end) = [];
% Xval_l(end) = [];
% Xdev{1} = Xtrain_l; Xdev{2} = Xval_l; Ydev{1} = Ytrain_l; Ydev{2} = Yval_l;
% Xdev{1} = toMat(Xtrain_l); Xdev{2} = toMat(Xval_l);
% Ydev{1}.Lfr = []; Ydev{2}.Lfr = [];
% for i = 1:20
%     Ydev{1}.Lfr = [Ydev{1}.Lfr Ytrain_l{i}.Lfr];
%     Ydev{2}.Lfr = [Ydev{2}.Lfr Yval_l{i}.Lfr];
% end

%% Baseline 
% First evaluation with euclidean distance
% profile -memory on
[~,S_eu,~] = g(params,Xdev{2},Ydev{2});
% profreport

%% Genetic DsTW algorithm 
% Evaluation function
if strcmp(params.scoreMeasure,'overlap')
    fEval = @(I) -fitnessFcn(I,X{1},Xdev{1},Ydev{1},Xdev{2},Ydev{2},params);    
elseif strcmp(params.scoreMeasure,'levenshtein')
    fEval = @(I) fitnessFcn(I,X{1},Xdev{1},Ydev{1},Xdev{2},Ydev{2},params);
end

% Display functions
fPlotComp = @(options,state,flag)plotMeanScores(options,state,flag,params,S_eu);
fPlotSI = @(options,state,flag)plotScoresPopul(options,state,flag,params);
fPlotSG = @(options,state,flag)plotScoreSegs(options,state,flag,params);

% Create function
fCreate = @(GenomeLength,FitnessFcn,options)createIniPopul(GenomeLength,FitnessFcn,options,X{1},params,seg0,fr_fixed);

% Mutation function
fMutation = @(parents,options,nvars,FitnessFcn,state,thisScore,...
    thisPopulation)mutationFcn(parents,options,nvars,FitnessFcn,...
    state,thisScore,thisPopulation,params,seg0,X{1});

% Crossover function
fCrossOver = @(parents,options,nvars,FitnessFcn,unused,thisPopulation)...
    crossOverFcn(parents,options,nvars,FitnessFcn,unused,thisPopulation,params);

% Options GA
lastGen = 4;
if exist(strcat('results/',DATATYPE,'/validation/Exp3/',params.Baseline,'Results',num2str(lastGen),'_',num2str(length(JOINTS)),COORDS,num2str(PERCENTDATA),'%_',num2str(NAT),'.mat'),'file')
    if strcmp(params.scoreMeasure,'overlap')
        load(strcat('results/',DATATYPE,'/validation/Exp3/',params.Baseline,'Results',num2str(lastGen),'_',num2str(length(JOINTS)),COORDS,num2str(PERCENTDATA),'%_',num2str(NAT),'.mat'));
    end   
    S = [];
    STATE = state;
    problem.nvars=size(state.Population,2);
    options = gaoptimset(options,'Generations',params.generations-lastGen);
    options = gaoptimset(options,'StallGenLimit',options.StallGenLimit+1);
else
    options = gaoptimset(@ga);
    options = gaoptimset(options,'PopulationSize',params.population);
    options = gaoptimset(options,'Generations',params.generations);
    options = gaoptimset(options,'TolCon',eps);
    options = gaoptimset(options,'StallGenLimit',params.generations);
    options = gaoptimset(options,'EliteCount',2);
    options = gaoptimset(options,'PlotFcns',{@gaplotscorediversity,@gaplotbestf,@gaplotrange,@gaplotdistance,fPlotComp,fPlotSI,fPlotSG});
    % options = gaoptimset(options,'PopInitRange',[params.k0 params.nmin*ones(1,100);params.N params.nmax*ones(1,100)]);
    options = gaoptimset(options,'CreationFcn',fCreate);
    options = gaoptimset(options,'MutationFcn',fMutation);
    options = gaoptimset(options,'CrossoverFcn',fCrossOver);
    options = gaoptimset(options,'Vectorized',params.vectorized);
    problem.nvars=1+params.N*2;
end
% Problem GA
problem.fitnessfcn=fEval;
problem.options=options;

% Run GA
tic;
if strcmp(params.vectorized,'on')
    if matlabpool('size') > 0
        matlabpool close force;
    end    
    pool(3);         % Cache with the populations
end
[x,finalOverlap,exitFlag,output,population,scores] = ga(problem);
if strcmp(params.vectorized,'on')
    matlabpool close;        % Cache with the populations
end
toc;
