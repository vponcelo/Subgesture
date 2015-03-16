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
COORDS = 'world'; NAT = 3;
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
params.phmm.Dtrain = cell(1,params.phmm.folds);
params.phmm.hmmTR_f = cell(1,params.phmm.folds); params.phmm.hmmE_f = cell(1,params.phmm.folds);
params.phmm.model = cell(1,params.phmm.folds);
params.phmm.pTrain_f = cell(1,params.phmm.folds); params.phmm.pVal_f = cell(1,params.phmm.folds);
params.phmm.C = cell(1,params.phmm.folds); params.phmm.SM = cell(1,params.phmm.folds);
params.phmm.map = cell(1,params.phmm.folds); params.phmm.minProb = cell(1,params.phmm.folds);
params.phmm.path = cell(1,params.phmm.folds);

%% Train and save learning results
if ~exist(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'),'file'),
    
    for k = 1:params.phmm.folds,
        display(sprintf('\n Fold %d',k));        
        params.phmm.Dtrain{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.hmmTR_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.hmmE_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.model{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.C{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.SM{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.pTrain_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.pVal_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.mapHMMtrain{k} = zeros(1,length(Xtrain_l)-1);
        params.phmm.mapHMMval{k} = zeros(1,length(Xtrain_l)-1);
        params.phmm.minProb{k} = zeros(1,length(Xtrain_l)-1);
        
        %% Training HMM models
        for l = 1:length(Xtrain_l)-1
            % Obtain Learning data
            Ytrain = Ydev{1}.Lfr(Ydev{1}.Lfr==l); Xtrain = Xtrain_l{l};     % gesture cells
%             Xtrain = Xdev{1}(Ydev{1}.Lfr==l,:);  
            if ~strcmp(params.phmm.clustType,'none')
                % Get data clusters (baseline 3)
                display('Discretizing training sequences in Key Poses ...');
                params.phmm.C{k}{l} = performClustering(Xtrain,Ytrain,params.phmm.clustType,params.phmm.kD,params.phmm.cIters);
                Xtrain = discretizeSequence(params.phmm.C{k}{l},Xtrain);
            end
            
            if strcmp(params.phmm.varType,'discrete')                
                % Get subgestures from training and validation
                if params.phmm.hmm
                    % Obtain Subgesture Model for training/learning data
                    [~,~,mErrsV,~,timeV,~,Z] = runKMeansDTW(params.version,params.k0,'dtwCost',params.k0,[],[],[],Ytrain,[],Xtrain,[]);
                    [~,kV] = min(mErrsV);
                    params.phmm.SM{k}{l} = Z{kV}{timeV}; emptyCells = cellfun(@isempty,params.phmm.SM{k}{l}); params.phmm.SM{k}{l}(emptyCells) = [];
                    display('Computing the costs of the training sequences in terms of SM and discretizing to the minimum cost ... ');
                    if iscell(Xtrain)
                        params.phmm.Dtrain{k}{l} = cell(1,length(Xtrain));
                        for sample = 1:length(Xtrain)
                            KM = getUpdatedCosts(Xtrain{sample},params.phmm.SM{k}{l});
                            [~,params.phmm.Dtrain{k}{l}{sample}] = min(KM);
%                             params.phmm.Dtrain{k}{l}{sample} = Xtrain{sample};
                        end
%                         params.phmm.Dtrain{k}{l} = cell2mat(params.phmm.Dtrain{k}{l});
                    else
                        KM = getUpdatedCosts(Xtrain,params.phmm.SM{k}{l});
                        [~,params.phmm.Dtrain{k}{l}] = min(KM);
%                         params.phmm.Dtrain{k}{l} = Xtrain;
                    end
                end
            end
            if strcmp(params.phmm.varType,'discrete')
                display(sprintf('Learning the HMM Model for gesture %d ...\n',l));                
                [params.phmm.hmmTR_f{k}{l}, params.phmm.hmmE_f{k}{l}] = ...
                    learnHMM(params.phmm.states,params.phmm.Dtrain{k}{l},params.phmm.it);
            elseif strcmp(params.phmm.clustType,'none')
                if strcmp(params.phmm.varType,'discrete')
                    [params.phmm.model{k}{l}, params.phmm.phmmloglikHist] = hmmFit(params.phmm.Dtrain{k}{l}, params.phmm.states, params.phmm.varType);
                elseif strcmp(params.phmm.varType,'gauss')
                    [params.phmm.model{k}{l}, params.phmm.phmmloglikHist] = hmmFit(Xtrain, params.phmm.states, params.phmm.varType);
                elseif strcmp(params.phmm.varType,'mixgausstied')
                    [~,kc] = kmeans(Xtrain,params.phmm.states);
                    [params.phmm.model{k}{l}, params.phmm.loglikHist] = hmmFit(Xtrain, params.phmm.states, params.phmm.varType, 'nmix', kc);
                end
            else
                error('testHMM:hmmTrainError','Error on the HMM training settings. Check varType and clustType parameters');
            end
        end
        
        %% Evaluating sequences for each model learnt
        for l = 1:length(Xtrain_l)-1            
            % Xval = Xdev{2}(Ydev{2}.Lfr==l,:);    % gesture vector
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
                    Xval = Xval(Yval == l,:);  % baseline 1+2: get current gesture label
                end
            else
                Yval = Ydev{2}.Lfr(Ydev{2}.Lfr==l); Xval = Xval_l{l};                   % gesture cells
            end
            if ~strcmp(params.phmm.clustType,'none')
                %% Get data clusters (baseline 3)
                display('Discretizing training sequences in Key Poses ...');
                Xval = discretizeSequence(params.phmm.C{k}{l},Xval);
            end
            
            % 1)    Testear todas las seqs de validación de un gesto con
            %       cada HMM aprendida. Tomar, por cada seq testeada, el
            %       likelihood y ver si pertenece a la HMM del gesto
            %       correspondiente.
            % 1Bis) Ver si son similares los resultados haciendo lo mismo testeando las secuencias de training            
            % 2)    Ver si la discretización por temporal clustering mejora
            % 3)    Ver cómo se tratan los estados de las secuencias y si 
            %       hay que añadir inicio-fin.
            % 4)    Comparar con la otra implementación.
            
            if strcmp(params.phmm.varType,'discrete')
                %% Get subgestures from training and validation
                if params.phmm.hmm
                    display('Computing the costs of the validation sequences in terms of SM and discretizing to the minimum cost ... ');
                    Dval = cell(1,length(Xval));
                    if iscell(Xval)
                        for sample = 1:length(Xval)
                            KT = getUpdatedCosts(Xval{sample},params.phmm.SM{k}{l});     
                            [~,Dval{sample}] = min(KT);
%                             Dval{sample} = Xval{sample};
                        end
                    else
                        KT = getUpdatedCosts(Xval,params.phmm.SM{k}{l});
                        [~,Dval] = min(KT);
%                         Dval = Xval;
                    end
                end
            end
            
            params.phmm.pTrain_f{k}{l} = zeros(length(params.phmm.Dtrain{k}),length(params.phmm.Dtrain{k}{l}));
            if iscell(Xval), 
                params.phmm.pVal_f{k}{l} = zeros(length(Xtrain_l)-1,length(Dval));
            else
                params.phmm.pVal_f{k}{l} = zeros(length(Xtrain_l)-1,1);
            end
            
            display(sprintf('Evaluation of training sequences of gesture %d with each HMM model ...',l));
            for m = 1:length(Xtrain_l)-1
                if strcmp(params.phmm.varType,'discrete')
                    if iscell(params.phmm.Dtrain{k}{l})
%                         display(sprintf('Evaluation of training sequences of gesture %d with the HMM model of gesture %d ...',l,m));
                        params.phmm.pTrain_f{k}{l}(m,:) = evaluateSequences([],params.phmm.Dtrain{k}{l},params.phmm.hmmTR_f{k}{m},params.phmm.hmmE_f{k}{m});
                    end
%                     display(sprintf('Evaluation of validation sequences of gesture %d with the HMM model of gesture %d ...',l,m));
                    params.phmm.pVal_f{k}{l}(m,:) = evaluateSequences([],Dval,params.phmm.hmmTR_f{k}{m},params.phmm.hmmE_f{k}{m});                    
                elseif strcmp(params.phmm.clustType,'none')
                    if strcmp(params.phmm.varType,'discrete')
                        params.phmm.path{k}{l} = hmmMap(params.phmm.model{k}{l}, Dval);
                    elseif strcmp(params.phmm.varType,'gauss')
                        params.phmm.path{k}{l} = hmmMap(params.phmm.model{k}{l}, Xval_l{l});
                    elseif strcmp(params.phmm.varType,'mixgausstied')
                        params.phmm.path{k}{l} = hmmMap(params.phmm.model{k}{l}, Xval_l{l});
                    end
                else
                    error('testHMM:hmmTrainError','Error on the HMM training settings. Check varType and clustType parameters');
                end
            end
            [~,params.phmm.mapHMMtrain{k}(l)] = max(mean(params.phmm.pTrain_f{k}{l}'));
            if size(params.phmm.pVal_f{k}{l},2) > 1
                [~,params.phmm.mapHMMval{k}(l)] = max(mean(params.phmm.pVal_f{k}{l}'));
            else
                [~,params.phmm.mapHMMval{k}(l)] = max(params.phmm.pVal_f{k}{l});
            end
            mapHMM = params.phmm.mapHMMval{k}(l);
            
            % Plot results of the model showing learning and predictive capabilities 
%                 plotResults(params.phmm.pTrain_f{k}{mapHMM},params.phmm.pVal_f{k}{mapHMM},...
%                     params.phmm.hmmE_f{k}{mapHMM},params.phmm.states,k);

            %% minimum probability to define the threshold
            params.phmm.minProb{k}(l) = min(params.phmm.pVal_f{k}{l}(mapHMM,:));
        end
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
