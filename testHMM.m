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
    Indices{l} = cell(1,params.phmm.folds);
    Indices{l} = crossvalind('Kfold',length(Xtrain_l{l}),params.phmm.folds);        
end
params.phmm.C = cell(1,params.phmm.folds);
params.phmm.hmmTR_f = cell(1,params.phmm.folds); params.phmm.hmmE_f = cell(1,params.phmm.folds);
params.phmm.pTrain_f = cell(1,params.phmm.folds); params.phmm.pVal_f = cell(1,params.phmm.folds);

%% Train and save learning results
if ~exist(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'),'file'),
    for k = 1:params.phmm.folds,
        display(sprintf('\n Fold %d',k));
        
        params.phmm.hmmTR_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.hmmE_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.pTrain_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.pVal_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.C{k} = cell(1,length(Xtrain_l)-1);
        for l = 1:length(Xtrain_l)-1

            %% Obtain Training data
            if params.phmm.folds > 1
                % grouped by gesture class
                [Xtrain,Ytrain] = getTrainingData(Xtrain_l,k,Indices,Ydev{1},l);                             
            else
                Xtrain = Xtrain_l{l}; Ytrain = Ydev{1}.Lfr(Ydev{1}.Lfr==l);                
            end            

            if strcmp(params.phmm.varType,'discrete')
                %% Get data clusters
                params.phmm.C{k}{l} = performClustering(Xtrain,Ytrain,params.phmm.clustType,params.phmm.kD,params.phmm.cIters);

                %% Discretize Training and Test data        
                Dtrain = discretizeData(params.phmm.C{k}{l},Xtrain);
                Dval = discretizeData(params.phmm.C{k}{l},Xdev{2});
            end
            if ~strcmp(params.phmm.clustType,'none') && strcmp(params.phmm.varType,'discrete')
                %% Test number of states 
                [params.phmm.hmmTR_f{k}{l},params.phmm.hmmE_f{k}{l},params.phmm.states,...
                    params.phmm.pTrain_f{k}{l},params.phmm.pVal_f{k}{l}] = ...
                    learnEvalModel(Dtrain,params.phmm.C{k}{l},Xtrain,Xval_l{l},params.phmm.it,params.phmm.states);
                
                %% Plot results of the model showing learning and predictive capabilities 
                plotResults(params.phmm.pTrain_f{k}{l},params.phmm.pVal_f{k}{l},...
                    params.phmm.hmmE_f{k}{l},params.phmm.states,k);

                %% minimum probability to define the threshold                
                %minModelProb{k}{l} = min(pVal);
            elseif strcmp(params.phmm.clustType,'none')
                if strcmp(params.phmm.varType,'discrete') 
                    tic;
                    [params.phmm.model, params.phmm.phmmloglikHist] = hmmFit(Xtrain, params.phmm.states, params.phmm.varType);
                    toc;
                    tic;
                    path = hmmMap(model, Dval);
                    toc;
                elseif strcmp(params.phmm.varType,'gauss')
                    tic;
                    [params.phmm.model, params.phmm.phmmloglikHist] = hmmFit(Xtrain, params.phmm.states, params.phmm.varType);
                    toc;
                    tic;
                    path = hmmMap(model, Xval_l{l});
                    toc;
                elseif strcmp(params.phmm.varType,'mixgausstied')
                    tic;
                    [~,kc] = kmeans(Xtrain,params.phmm.states);
                    [params.phmm.model, params.phmm.loglikHist] = hmmFit(Xtrain, params.phmm.states, params.phmm.varType, 'nmix', kc);                    
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
    save(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'),'params');
    display('Done!');
else
    display('Showing Learning results for each fold ...');
    load(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'));
    minModelProb = cell(1,params.phmm.folds);
    minP = ones*inf(1,params.phmm.folds); posG = ones*inf(1,params.phmm.folds); 
    for k = 1:params.phmm.folds,
        minModelProb{k} = zeros(1,length(Xtrain_l)-1);
        for l = 1:length(Xtrain_l)-1
            minModelProb = cell(1,params.phmm.folds);
            plotResults(params.phmm.pTrain_f{k}{l},params.phmm.pVal_f{k}{l},...
                params.phmm.hmmE_f{k}{l},params.phmm.states,k);
            minModelProb{k}(l) = min(params.phmm.pVal_f{k}{l});
        end
        [minP(k),posG(k)] = min(minModelProb{k});
    end
    display('Done!');

    %% Evaluate Test data
    [threshold,k] = min(minP);
    if threshold < 0.5
        threshold = 0.5; 
    elseif threshold > 0.8
        threshold = 0.8;
    end
    minModelProbs = params.phmm.pVal_f{k}{posG(k)};

    display('Evaluating the final Model with the test set...');
    if ~exist(strcat('results/',DATATYPE,'/validation/hmm/resProbTestSeqs.mat'),'file'),
        testProbs=evaluateSequences(params.phmm.C{k}{posG(k)},Xtest,...
            params.phmm.hmmTR_f{k}{posG(k)},params.phmm.hmmE_f{k}{posG(k)});
        plotResults(-1,testProbs,params.phmm.hmmE_f{k}{posG(k)},params.phmm.states,k);

        hits = sum(testProbs > threshold);
        accuracy = hits/length(testProbs);

        save(strcat('results/',DATATYPE,'/validation/hmm/resProbTestSeqs.mat'),'testProbs','minModelProbs','threshold','accuracy');
        display('Done!');
    else
        load(strcat('results/',DATATYPE,'/validation/hmm/resProbTestSeqs.mat'));
        plotResults(-1,testProbs,params.phmm.hmmE_f{k},params.phmm.states,k);
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



