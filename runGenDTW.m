function runGenDTW(measure,lastGen)

close all
% clear all
addPath

%% Generate global variables for the GA, cache and parameters
varload
if nargin == 0
    measure = 'overlap';
    lastGen = 0;
end
switch params.score2optim
    case 'o', if ~params.classification, params.score2optim = 1; else error('g:optErr','Spotting is not allowed in classification ...for now... (ODL!)'); end
    case 'p', if ~params.classification, params.score2optim = 2; else params.score2optim = 1; end
    case 'r', if ~params.classification, params.score2optim = 3; else params.score2optim = 2; end
    case 'a', if ~params.classification, params.score2optim = 4; else params.score2optim = 3; end
end
params.scoreMeasure = measure;  % Score Measure: 'overlap' or 'levenshtein'
if strcmp(measure,'overlap')
    CACHE.eval = zeros(...      % Cache with the evaluations of each individual
        params.population*100,1);
    if params.classification, nsc = 3; else nsc = 4; end
    CACHE.scores = zeros(...      % Cache with the scores of each individual
        params.population*100,nsc);
elseif strcmp(measure,'levenshtein')
    CACHE.eval = inf*ones(...   % Cache with the evaluations of each individual
        params.population*100,1);
    CACHE.eval = inf*ones(...   % Cache with the scores of each individual
        params.population*100,4);
end

%% Prepare training data depending on the chosen option and parameters
% Load data:
%     if nframesSeg is 0, then initial segmentation is generated from the skeleton labels
% [X,Y,Xtest,Ytest] = prepareData(nrsamples,nseqs,nframesSeg,params.k0);
% display('Press a key to continue...');
% pause();

%% Compute initial segmentation from motion
seg0 = [];
% [seg0,fr_fixed,params] = computeMagnitudes(X{1},params);

%% Prepare training data depending on the chosen option and parameters
NORMTYPE = 'none'; COORDS = 'world'; NAT = 3;
% Load data:
%     if nframesSeg is 0, then initial segmentation is generated from the skeleton labels
[X,Y,Xtest,Ytest] = prepareData(nrsamples,nseqs,nframesSeg,params.k0);
% X = setDerivatives(X);
% display('Press a key to continue...');
% pause();

%% Obtain all samples grouped (labeled) by gestures
%Xtrain_l = getGroupedGestures(X,Y,0);
Xtrain_l = getGroupedGestures(X,Y,1,[],[]); if sum(cellfun(@isempty,Xtrain_l)), error('Empty gesture classes'); end
Xval_l = getGroupedGestures(X,Y,2,[],[]); if sum(cellfun(@isempty,Xval_l)), error('Empty gesture classes'); end
if ~isempty(Xtest) && ~isempty(Ytest)
    Xtest_l = getGroupedGestures(X,Y,3,Xtest,Ytest); if sum(cellfun(@isempty,Xtest_l)), error('Empty gesture classes'); end
end

%% Compute median models from training/learning data
% profile -memory on
if strcmp(DATATYPE,'msr3dS') || strcmp(DATATYPE,'msr3d1') || strcmp(DATATYPE,'msr3d2') || ...
    strcmp(DATATYPE,'msr3d3') || strcmp(DATATYPE,'msr3d4') || strcmp(DATATYPE,'msr3d5') || strcmp(DATATYPE,'msract3d'), 
    nModels = length(Xtrain_l); 
else nModels = length(Xtrain_l)-1; 
end  % -1 indicates don't consider iddle gesture
if ~params.phmm.hmm, [params.M,params.lmodel] = getModels(Xtrain_l,nModels,params); end
% profreport

%% Generate development sequences
% l = [24 78 150];    % 78 (more samples for each gesture when k=3);
l = [];
[Xdev,Ydev,Xtest,Ytest] = getDevSequences(X,Y,Xtest,Ytest,l,noise,secsBatch,nSampGest);
clear X Y

%% Baseline 
% First evaluation with euclidean distance
% profile -memory on
S_base = 0;
if params.darwin
%     S_base = testDarwin(Xdev{1},Xdev{2},Ydev{1},Ydev{2});
    S_base = testDarwin(Xdev{1},Xtest,Ydev{1},Ytest);
else
    if ~params.phmm.hmm
        if params.classification
            [model,S_base,bestScores,~] = g(params,Xval_l,Ydev{2});
            [~,S_base,bestScores,~] = g(model,Xtest_l,Ytest);
        else
            [model,S_base,bestScores,~] = g(params,Xdev{2},Ydev{2});
            [~,S_base,bestScores,~] = g(model,Xtest,Ytest);
        end
    else
        [S_base,model,bestScores] = testHMM(params);
        if params.classification
            Dtest = cell(1,length(Xtest_l));
            for l = 1:length(Xtest_l)
                Dtest{l} = cell(1,length(Xtest_l{l}));
                for sample = 1:length(Xtest_l{l})
                    KT = getUpdatedCosts(Xtest_l{l}{sample},model.SM);
                    [~,Dtest{l}{sample}] = min(KT);
                end
            end
        else
            KT = getUpdatedCosts(Xtest,model.SM);
            [~,Dtest] = min(KT);
        end
        [~,S_base,bestScores] = evalswHMM(model, Dtest, Ytest);
    end
end
S_base
% profreport

%% Genetic algorithm optimization
% Evaluation function
if strcmp(params.scoreMeasure,'overlap')
    fEval = @(I) -fitnessFcn(I,Xdev{1},Xtrain_l(1:nModels),Ydev{1},Xval_l(1:nModels),Xdev{2},Ydev{2},Xtest,Xtest_l,Ytest,params);    
elseif strcmp(params.scoreMeasure,'levenshtein') || params.phmm.hmm
    fEval = @(I) fitnessFcn(I,Xdev{1},Xtrain_l(1:nModels),Ydev{1},Xval_l(1:nModels),Xdev{2},Ydev{2},Xtest,Xtest_l,Ytest,params);
end

% Display functions
fPlotComp = @(options,state,flag)plotMeanScores(options,state,flag,params,S_base,Xtest,Ytest,Xtest_l);
fPlotSI = @(options,state,flag)plotScoresPopul(options,state,flag,params);
fPlotSG = @(options,state,flag)plotScoreSegs(options,state,flag,params);

% Create function
fCreate = @(GenomeLength,FitnessFcn,options)createIniPopul(GenomeLength,FitnessFcn,options,Xdev{1},params);

% Mutation function
fMutation = @(parents,options,nvars,FitnessFcn,state,thisScore,...
    thisPopulation)mutationFcn(parents,options,nvars,FitnessFcn,...
    state,thisScore,thisPopulation,params,seg0,Xdev{1});

% Crossover function
fCrossOver = @(parents,options,nvars,FitnessFcn,unused,thisPopulation)...
    crossOverFcn(parents,options,nvars,FitnessFcn,unused,thisPopulation,params,Xdev{1});

% Options GA
lastGen = 10;
if exist(strcat('results/',DATATYPE,'/validation/Exp3/gen',num2str(params.generations),'popul',num2str(params.population),'/',...
        params.Baseline,'_',params.mType,'_',num2str(lastGen),'gens','_',...
        num2str(length(JOINTS)),'joints',COORDS,'_','mod',num2str(NAT),'.mat'),'file')
    if strcmp(params.scoreMeasure,'overlap')
        load(strcat('results/',DATATYPE,'/validation/Exp3/gen',num2str(params.generations),'popul',num2str(params.population),'/',...
            params.Baseline,'_',params.mType,'_',num2str(lastGen),'gens','_',...
            num2str(length(JOINTS)),'joints',COORDS,'_','mod',num2str(NAT),'.mat'));
    end
    if ~isempty(STATE)
        state = STATE;
    elseif ~isempty(state)
        STATE = state;
    else
        error('runGenDTW:stateErr','State was not stored in previous execution');
    end
    if isempty(BESTIND)
        error('runGenDTW:bestIndErr','BestIndividual was not stored in previous execution');
    end
    bestind.model = []; bestind.state = []; BESTIND = [BESTIND bestind];
    problem.nvars=size(state.Population,2);
    options = gaoptimset(options,'PopulationSize',params.population);
    options = gaoptimset(options,'Vectorized',params.vectorized);
    options = gaoptimset(options,'MutationFcn',fMutation);
    options = gaoptimset(options,'CrossoverFcn',fCrossOver);
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
    if strcmp(params.mType,MEDIANTYPE{2}) || ...
            strcmp(params.mType,MEDIANTYPE{4}), problem.nvars = problem.nvars+2; 
    end
    
    % Cache with the populations
    CACHE.ind = int32(zeros(params.population*100,problem.nvars,'int32'));

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
    pool(2);         % Cache with the populations
end
[x,finalOverlap,exitFlag,output,population,scores] = ga(problem);
if strcmp(params.vectorized,'on')
    matlabpool close;        % Cache with the populations
end
toc;
