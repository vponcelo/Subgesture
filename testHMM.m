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

%% Obtain all training samples grouped (labeled) by gestures
Xtrain_l = getGroupedGestures(X,Y,1);
Xval_l = getGroupedGestures(X,Y,2);

%% Generate learning sequences
% l = [24 78 150];    % 78 (more samples for each gesture when k=3);
l = [];
[Xdev,Ydev] = getDevSequences(X,Y,l,noise,secsBatch,nSampGest);

%% Obtain Cross Validation subsets over training data
Indices = cell(1,length(Xtrain_l)-1);
for l = 1:length(Xtrain_l)-1
    Indices{l} = cell(1,phmm.folds);
    Indices{l} = crossvalind('Kfold',length(Xtrain_l{l}),phmm.folds);        
end
phmm.hmmTR_f = cell(1,phmm.folds); phmm.hmmE_f = cell(1,phmm.folds);
pTrain_f = cell(1,phmm.folds); pVal_f = cell(1,phmm.folds);

%% Train and save learning results
if ~exist(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'),'file'),
    for k = 1:phmm.folds,
        display(sprintf('\n Fold %d',k));
        
        phmm.hmmTR_f{k} = cell(1,length(Xtrain_l)-1);
        phmm.hmmE_f{k} = cell(1,length(Xtrain_l)-1);
        pTrain_f{k} = cell(1,length(Xtrain_l)-1);
        pVal_f{k} = cell(1,length(Xtrain_l)-1);
        for l = 1:length(Xtrain_l)-1

            %% Obtain Training data
            if phmm.folds > 1
                % grouped by gesture class
                [Xtrain,Ytrain] = getTrainingData(Xtrain_l,k,Indices,Ydev{1},l);                             
            else
                Xtrain = Xtrain_l{l}; Ytrain = Ydev{1}.Lfr(Ydev{1}.Lfr==l);                
            end            

            if strcmp(phmm.varType,'discrete')
                %% Get data clusters
                Ctrain = performClustering(Xtrain,Ytrain,phmm.clustType,phmm.kD,phmm.cIters);

                %% Discretize Training and Test data        
                Dtrain = discretizeData(Ctrain,Xtrain);
                Dval = discretizeData(Ctrain,Xdev{2});
            end
            if ~strcmp(phmm.clustType,'none') && strcmp(phmm.varType,'discrete')
                %% Test number of states 
                [phmm.hmmTR_f{k}{l},phmm.hmmE_f{k}{l},phmm.hmmStates,...
                    phmm.pTrain_f{k}{l},phmm.pVal_f{k}{l}] = ...
                    learnEvalModel(Dtrain,Ctrain,Xtrain,Xval_l{l},phmm.it,phmm.states);
                
                %% Plot results of the model showing learning and predictive capabilities 
                plotResults(phmm.pTrain_f{k}{l},phmm.pVal_f{k}{l},...
                    phmm.hmmE_f{k}{l},phmm.hmmStates,k);

                %% minimum probability to define the threshold                
                %minModelProb{k}{l} = min(pVal);
            elseif strcmp(phmm.clustType,'none')
                if strcmp(phmm.varType,'discrete') 
                    tic;
                    [phmm.model, phmm.phmmloglikHist] = hmmFit(Xtrain, phmm.states, phmm.varType);
                    toc;
                    tic;
                    path = hmmMap(model, Dval);
                    toc;
                elseif strcmp(phmm.varType,'gauss')
                    tic;
                    [phmm.model, phmm.phmmloglikHist] = hmmFit(Xtrain, phmm.states, phmm.varType);
                    toc;
                    tic;
                    path = hmmMap(model, Xval_l{l});
                    toc;
                elseif strcmp(phmm.varType,'mixgausstied')
                    tic;
                    [~,kc] = kmeans(Xtrain,phmm.states);
                    [phmm.model, phmm.loglikHist] = hmmFit(Xtrain, phmm.states, phmm.varType, 'nmix', kc);                    
                    toc;
                    tic;
                    path = hmmMap(model, Xval_l{l});
                    toc;
                end
                %[~,S_eu,~] = g(params,Xdev{2},Ydev{2});
                %S_eu
            else
                error('testHMM:hmmTrainError','Error on the HMM training settings. Check varType and clustType parameters');
            end
        end
    end
    display(sprintf('Saving model and results ...'));
    save(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'),'phmm');
    display('Done!');
else
    display('Showing Learning results for each fold ...');
    load(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'));
    for k = 1:params.folds,
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



