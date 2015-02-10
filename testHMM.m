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
% NAT = 3;
% NORMTYPE = 'none';
% COORDS = 'world';

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

%% Obtain all training samples grouped (labeled) by gestures
if ~phmm.wholeSeq
    Xtrain_l = getGroupedGestures(Xdev,Ydev,1);
else
    Xtrain_l = 1;
end

%% Obtain Cross Validation subsets over training data
if phmm.wholeSeq
    Indices = crossvalind('Kfold',length(Xdev{1}),phmm.folds);    
    phmm.hmmTR_f = cell(1,phmm.folds);
    phmm.hmmE_f = cell(1,phmm.folds);
    phmm.hmmStates = cell(1,phmm.folds);
    pVal_f = zeros(phmm.folds,round(length(X)/phmm.folds));
    pTrain_f = zeros(phmm.folds,round(size(pVal_f,2)*(phmm.folds-1)));
    minModelProb = ones(1,phmm.folds);
else
    Indices = cell(1,length(Xtrain_l));
    for l = 1:length(Xtrain_l)
        Indices{l} = cell(1,phmm.folds);
        Indices{l} = crossvalind('Kfold',length(Xtrain_l{l}),phmm.folds);        
    end    
end

%% Train and save learning results
if ~exist(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'),'file'),
    for k = 1:phmm.folds,        
        for l = 1:length(Xtrain_l)
            display(sprintf('\n Fold %d',k));

            %% Obtain Training data
            if phmm.folds > 1
                if ~iscell(Xtrain_l)
                    [Xtrain,Ytrain] = getTrainingData(Xdev{1},k,Indices,Ydev{1});
                else
                    % grouped by gesture class
                    [Xtrain,Ytrain] = getTrainingData(Xtrain_l,k,Indices,Ydev{1},l);             
                end
            else
                if ~iscell(Xtrain_l)
                    Xtrain = Xdev{1}; Ytrain = Ydev{1};                     
                else
                    Xtrain = Xtrain_l{l}; Ytrain = Ydev{1}.Lfr(Ydev{1}.Lfr==l);
                end
            end            

            if strcmp(phmm.varType,'discrete') && ~strcmp(clustType,'none')
                %% Get data clusters
                Ctrain = performClustering(Xtrain,Ytrain,phmm.clustType,phmm.kD,phmm.cIters);

                %% Discretize Training and Test data        
                Dtrain = discretizeData(Ctrain,Xtrain);
            
                %% Test number of states 
                display(sprintf('Learning the Model ...'));
                [hmmTR,hmmE,phmm.hmmStates,pTrain,pVal] = learnModel(Dtrain,Ctrain,Xtrain,Xdev{2},phmm.it,phmm.states);
                display(sprintf('Model learnt'));

                %% Plot results of the model showing learning and predictive capabilities 
                plotResults(pTrain,pVal,hmmE,hmmStates,k);

                %% save HMM transition and emissions, probability values for traning and validation data, and minimum probability to define the threshold
                phmm.hmmTR_f{k} = hmmTR;
                phmm.hmmE_f{k} = hmmE;
                phmm.hmmStates_f(k,:) = hmmStates;
                phmm.pTrain_f(k,:) = pTrain;
                phmm.pVal_f(k,:) = pVal;        
                minModelProb(k) = min(pVal);
            elseif strcmp(phmm.varType,'gauss') || strcmp(phmm.varType,'mixgausstied') || ...
                    strcmp(phmm.varType,'discrete') && strcmp(clustType,'none')
                
                if strcmp(phmm.varType,'mixgausstied')
                    tic;
                    [~,kMeansClusters] = kmeans(Xtrain,phmm.states);
                    [model, loglikHist] = hmmFit(Xtrain, phmm.states, phmm.varType, 'nmix', kMeansClusters);
                    toc;
                else
                    tic;
                    [model, loglikHist] = hmmFit(Xtrain, phmm.states, phmm.varType);
                    toc;
                end
                tic;
                path = hmmMap(model, Xdev{2});                
                toc;
                %[~,S_eu,~] = g(params,Xdev{2},Ydev{2});
                %S_eu
            else
                error('testHMM:hmmTrainError','Error on the HMM training settings. Check varType and clustType parameters');
            end
        end
    end
    display(sprintf('Saving results ...'));
    save(strcat('results/',DATATYPE,'/HMM/,learningResults.mat'),'phmm','pTrain_f','pVal_f','minModelProb');
    display('Done!');
else
    display('Showing Learning results for each fold ...');
    load('data/learningResults.mat');
    for k = 1:folds,
        plotResults(pTrain_f(k,:),pVal_f(k,:),phmm.hmmE_f{k},phmm.hmmStates,k);
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

        save('data/resProbTestSeqs.mat','testProbs','minModelProbs','threshold','accuracy','clustType');
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



