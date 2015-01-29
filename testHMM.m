function testHMM()

clear all
close all
addPath

%% Generate global variables, cache and parameters
varload
if nargin == 0
    measure = 'overlap';
    lastGen = 0;
end
params.scoreMeasure = measure;  % Score Measure: 'overlap' or 'levenshtein'
clear CACHE S BASELINE OPTIONS STATE;


%% Prepare training data depending on the chosen option and parameters
% Load data:
%     if nframesSeg is 0, then initial segmentation is generated from the skeleton labels
[X,Y,Xtest,Ytest] = prepareData(nrsamples,nseqs,nframesSeg,params.k0);
% display('Press a key to continue...');
% pause();

%% Generate learning sequences
% l = [24 78 150];    % 78 (more samples for each gesture when k=3);
l = [];
[Xdev,Ydev] = getDevSequences(X,Y,l,noise,secsBatch,nSampGest);

%% Obtain Cross Validation subsets over training data
if folds > 1
    Indices = crossvalind('Kfold',length(Xdev{1}),folds);
else
    Xtrain = Xdev{1}; Ytrain = Ydev{1}.Lfr;
end
hmmTR_f = cell(1,folds);
hmmE_f = cell(1,folds);
pVal_f = zeros(folds,round(length(X)/folds));
pTrain_f = zeros(folds,round(size(pVal_f,2)*(folds-1)));
minModelProb = ones(1,folds);

%% Train and save learning results
if ~exist(strcat('results/',DATATYPE,'/HMM/,learningResults.mat'),'file'),
    for k = 1:folds,

        display(sprintf('\n Fold %d',k));

        %% Obtain Training data
        if folds > 1
            [Xtrain,Ytrain] = getTrainingData(Xdev{1},k,Indices,Ydev{1}); % supervised 
        end

        %% Get data clusters
        Ctrain = performClustering(Xtrain,Ytrain,clustType,numClusters,numIterations);

        %% Discretize Training and Test data        
        Dtrain = discretizeData(Ctrain,Xtrain);
       
        %% Test number of states 
        display(sprintf('Learning the Model ...'));        
        [hmmTR,hmmE,hmmStates,pTrain,pVal] = learnModel(Dtrain,Ctrain,Xtrain,Xdev{2},hmmIters);

        %% Plot results of the model showing learning and predictive capabilities 
        plotResults(pTrain,pVal,hmmE,hmmStates,k);

        %% save HMM transition and emissions, probability values for traning and validation data, and minimum probability to define the threshold
        hmmTR_f{k} = hmmTR;
        hmmE_f{k} = hmmE;
        pTrain_f(k,:) = pTrain;
        pVal_f(k,:) = pVal;
        minModelProb(k) = min(pVal);
    end
    save(strcat('results/',DATATYPE,'/HMM/,learningResults.mat'),'hmmE_f','pTrain_f','pVal_f','folds','hmmStates','minModelProb','Ctrain','datatype','clustType');
    display('Done!');
else
    display('Showing Learning results for each fold ...');
    load('data/learningResults.mat');
    for k = 1:folds,
        plotResults(pTrain_f(k,:),pVal_f(k,:),hmmE_f{k},hmmStates,k);
    end
    display('Done!');

    %% Evaluate Test data

    [threshold,k] = min(minModelProb);
    if threshold < 0.5
        threshold = 0.5; 
    elseif threshold > 0.8
        threshold = 0.8; 
    end
    minModelProbs = pVal_f(k,:);

    display('Evaluating the final Model with the test set...');
    if ~exist('data/resProbTestSeqs.mat','file'),
        testProbs=evaluateSequences(Ctrain,Xtest,hmmTR_f{k},hmmE_f{k});
        plotResults(-1,testProbs,hmmE_f{k},hmmStates,k);

        hits = sum(testProbs > threshold);
        accuracy = hits/length(testProbs);

        save('data/resProbTestSeqs.mat','testProbs','minModelProbs','threshold','accuracy','datatype','clustType');
        display('Done!');
    else
        load('data/resProbTestSeqs.mat');
        plotResults(-1,testProbs,hmmE_f{k},hmmStates,k);
    end

    figure,
    hold on
    plot(minModelProbs);
    title('Blue: validation set of the selected fold. Red: Test set');
    plot(testProbs,'red');
    ylabel('Probability values');
    xlabel('Sequence number');
    hold off

    display('Done!');
end



