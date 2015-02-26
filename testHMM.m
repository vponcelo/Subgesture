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
params.phmm.map = cell(1,params.phmm.folds); params.phmm.minProb = cell(1,params.phmm.folds);

%% Train and save learning results
if ~exist(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'),'file'),
    for k = 1:params.phmm.folds,
        display(sprintf('\n Fold %d',k));        
        params.phmm.hmmTR_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.hmmE_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.pTrain_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.pVal_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.C{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.map{k} = zeros(1,length(Xtrain_l)-1);
        for l = 1:length(Xtrain_l)-1
            %% Obtain Learning data
            Ytrain = Ydev{1}.Lfr(Ydev{1}.Lfr==l); Xtrain = Xtrain_l{l};                 % gesture cells
%             Xtrain = Xdev{1}(Ydev{1}.Lfr==l,:);  Xval = Xdev{2}(Ydev{2}.Lfr==l,:);    % gesture vector
            if params.sw > 0
                Xval = [];
                while isempty(Xval)     % baseline 2 | 1+2
                    r = inf;
                    if params.sw == length(Ydev{2}.Lfr)
                        r = 1;
                    else
                        while any(r > length(Ydev{2}.Lfr)-params.sw)
                            r=randperm(round(length(Ydev{2}.Lfr)),1);
                        end
                    end
                    seg=r(1):min(r(1)+params.sw,length(Ydev{2}.Lfr));
                    Xval=Xdev{2}(seg,:); Yval=Ydev{2}.Lfr(seg);                         % gesture vector
%                     Xval = Xval(Yval == l,:);  % baseline 1+2: get current gesture label
                end
            else
                Yval = Ydev{2}.Lfr(Ydev{2}.Lfr==l); Xval = Xval_l{l};                   % gesture cells
            end

            if strcmp(params.phmm.varType,'discrete')
                if ~strcmp(params.phmm.clustType,'none')
                    %% Get data clusters (baseline 3)
                    % s'ha de representar X en bokps fent clustering
                    Ctrain = performClustering(Xtrain,Ytrain,params.phmm.clustType,params.phmm.kD,params.phmm.cIters);
                    Xtrain = discretizeSequence(Ctrain,Xtrain);
                    Xval = discretizeSequence(Ctrain,Xval);
                end
                %% Get subgestures from training and validation
                if params.phmm.hmm
                    %% Obtain Subgesture Model for training/learning data
                    [~,~,mErrsV,~,timeV,~,Z] = runKMeansDTW(params.version,params.k0,'dtwCost',params.k0,[],[],[],Ytrain,[],Xtrain,[]);
                    [~,kV] = min(mErrsV);
                    SM = Z{kV}{timeV};
                    emptyCells = cellfun(@isempty,SM);
                    SM(emptyCells) = [];
%                     display('Computing the costs of the training sequences in terms of SM and discretizing to the minimum ... ');
                    if iscell(Xtrain)
                        Dtrain = cell(1,length(Xtrain));
                        for sample = 1:length(Xtrain)
                            KM = getUpdatedCosts(Xtrain{sample},SM);
                            [~,Dtrain{sample}] = min(KM);
                        end
                        Dtrain = cell2mat(Dtrain);
                    else
                        KM = getUpdatedCosts(Xtrain,SM);
                        [~,Dtrain] = min(KM);
                    end
%                     display('Computing the costs of the validation sequences in terms of SM and discretizing to the minimum ... ');
                    Dval = cell(1,length(Xval));
                    if iscell(Xval)
                        for sample = 1:length(Xval)
                            KT = getUpdatedCosts(Xval{sample},SM);                            
                            [~,Dval{sample}] = min(KT);
                        end
                    else
                        KT = getUpdatedCosts(Xval,SM);
                        [~,Dval] = min(KT);
                    end
                end
            end
            if strcmp(params.phmm.varType,'discrete')
                display(sprintf('Learning the and evaluating the HMM Model for gesture %d ...',l));                
                [params.phmm.hmmTR_f{k}{l}, params.phmm.hmmE_f{k}{l}] = ...
                    learnHMM(params.phmm.states,Dtrain,params.phmm.it);
                if iscell(Dtrain)
%                     display('Evaluation of training ...');
                    params.phmm.pTrain_f{k}{l} = evaluateSequences([],Dtrain,params.phmm.hmmTR_f{k}{l},params.phmm.hmmE_f{k}{l});
%                     max(params.phmm.pTrain_f{k}{l})
                end
%                 display('Evaluation of validation ...');
                params.phmm.pVal_f{k}{l} = evaluateSequences([],Dval,params.phmm.hmmTR_f{k}{l},params.phmm.hmmE_f{k}{l});
                params.phmm.map{k}(l) = max(params.phmm.pVal_f{k}{l});
%                 params.phmm.map{k}(l)
                
                %% Plot results of the model showing learning and predictive capabilities 
%                 plotResults(params.phmm.pTrain_f{k}{l},params.phmm.pVal_f{k}{l},...
%                     params.phmm.hmmE_f{k}{l},params.phmm.states,k);

                %% minimum probability to define the threshold
                params.phmm.minProb{k} = min(params.phmm.map{k}(l));
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
            else
                error('testHMM:hmmTrainError','Error on the HMM training settings. Check varType and clustType parameters');
            end
        end
        mean(params.phmm.map{k})
    end
    display(sprintf('Saving model and results ...'));
    save(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'),'params');
    display('Done!');
else
    display('Showing Learning results for each fold ...');
    load(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'));
    for k = 1:params.phmm.folds,
        for l = 1:length(Xtrain_l)-1
            plotResults(params.phmm.pTrain_f{k}{l},params.phmm.pVal_f{k}{l},...
                params.phmm.hmmE_f{k}{l},params.phmm.states,k);            
        end
        [minP(k),posG(k)] = min(params.phmm.minProb{k});
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
