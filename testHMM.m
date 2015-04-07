function [score,model] = testHMM(params)

% clear all
% close all
% addPath

%% Generate global variables, cache and parameters
varload
if ~exist('params','var')
    params.scoreMeasure = 'overlap';  % Score Measure: 'overlap' or 'levenshtein'    
    clear CACHE S BASELINE OPTIONS STATE;
else
    nSampGest = 10;
end

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
Xval = Xdev{2}; Yval = Ydev{2};

%% Obtain Cross Validation subsets over training data
Indices = cell(1,length(Xtrain_l)-1);
for l = 1:length(Xtrain_l)-1
    Indices{l} = cell(1,params.phmm.folds);
    Indices{l} = crossvalind('Kfold',length(Xtrain_l{l}),params.phmm.folds);        
end
if ~strcmp(params.phmm.clustType,'none') 
    params.phmm.C = cell(1,params.phmm.folds);
end
params.phmm.SM = cell(1,params.phmm.folds);
if strcmp(params.phmm.varType,'discrete') 
    params.phmm.Dtrain = cell(1,params.phmm.folds);
end
params.phmm.hmmTR_f = cell(1,params.phmm.folds); params.phmm.hmmE_f = cell(1,params.phmm.folds);
params.phmm.model = cell(1,params.phmm.folds); params.phmm.path = cell(1,params.phmm.folds);
params.phmm.pTrain_f = cell(1,params.phmm.folds); params.phmm.pVal_f = cell(1,params.phmm.folds); 
params.phmm.mapHMMtrain = cell(1,params.phmm.folds); params.phmm.mapHMMval = cell(1,params.phmm.folds); params.phmm.minProb = cell(1,params.phmm.folds); 
params.phmm.accTrain = cell(1,params.phmm.folds); params.phmm.accLearn = cell(1,params.phmm.folds);

%% Train and save learning results
if ~exist(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'),'file'),
    
    for k = 1:params.phmm.folds,
        display(sprintf('\n Fold %d',k));
        if strcmp(params.phmm.varType,'discrete') 
            params.phmm.Dtrain{k} = cell(1,length(Xtrain_l)-1);
        end
        params.phmm.pTrain_f{k} = cell(1,length(Xtrain_l)-1); params.phmm.pVal_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.mapHMMtrain{k} = cell(1,length(Xtrain_l)-1); params.phmm.mapHMMval{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.minProb{k} = zeros(1,length(Xtrain_l)-1);
        params.phmm.accTrain{k} = zeros(1,length(Xtrain_l)-1); params.phmm.accLearn{k} = zeros(1,length(Xtrain_l)-1);
        params.phmm.hmmTR_f{k} = cell(1,length(Xtrain_l)-1); params.phmm.hmmE_f{k} = cell(1,length(Xtrain_l)-1);
        params.phmm.model{k} = cell(1,length(Xtrain_l)-1); params.phmm.path{k} = cell(1,length(Xtrain_l)-1);
        
        %% Select a sample containing gesture samples of each class
        Xtrain = [];
        for l = 1:length(Xtrain_l)
            if ~nSampGest
                Xtrain = [Xtrain Xtrain_l{l}];
            elseif nSampGest > 0
                r = randi([1,length(Xtrain_l{l})-nSampGest]);
                Xtrain = [Xtrain Xtrain_l{l}(r:r+(nSampGest-1))];
            end            
        end
        params.phmm.kD = length(Xtrain);
        if ~strcmp(params.phmm.clustType,'none')
            %% Clustering over the sample
            params.phmm.C{k} = ...
            performClustering(Xtrain,[],params.phmm.clustType,params.phmm.kD,params.phmm.cIters);
        end
        if strcmp(params.phmm.varType,'discrete')   % comment this case to avoid temporal clustering
            %% Temporal Clustering over the sample
            Xt = [];
            for gest = 1:length(Xtrain)
                Xsplit = cell(1,params.phmm.states);
                div = round(size(Xtrain{gest},1)/params.phmm.states);
                s = 1;
                for i = 1:params.phmm.states
                    e = s+div;
                    Xsplit{i} = Xtrain{gest}(s:min(e,size(Xtrain{gest},1)),:);
                    s = e+1;
                end
                Xt = [Xt Xsplit];
            end
            emptyCells = cellfun(@isempty,Xt); Xt(emptyCells) = [];
            k0 = round(params.phmm.kD/2);
            [~,~,mErrsV,~,timeV,~,Z] = runKMeansDTW(params,k0,k0,[],[],[],[],[],Xt,[]);
            [~,kV] = min(mErrsV);
            params.phmm.SM{k} = Z{kV}{timeV}; emptyCells = cellfun(@isempty,params.phmm.SM{k}); params.phmm.SM{k}(emptyCells) = [];
        end
        
        % Training HMM models
        for l = 1:length(Xtrain_l)-1
            %% Obtain Training data for each gesture class            
            Xtrain = Xtrain_l{l}(1:1+(nSampGest-1));
%             Ytrain = Ydev{1}.Lfr(Ydev{1}.Lfr==l); Xtrain = Xtrain_l{l}; % gesture cells
%             Xtrain = Xdev{1}(Ydev{1}.Lfr==l,:);                         % gesture vector
            %% Discretize training data
            if strcmp(params.phmm.varType,'discrete')
                if ~strcmp(params.phmm.clustType,'none')
                    display('Discretizing training sequence/s in Key Poses ...');
                    Xtrain = discretizeSequence(params.phmm.C{k},Xtrain);
                end
                if iscell(Xtrain)
                    params.phmm.Dtrain{k}{l} = cell(1,length(Xtrain));
                    for sample = 1:length(Xtrain)
                        if ~isempty(params.phmm.SM{k})
                            %% Obtain Subgesture cost representations for training data
                            if sample==1,display('Computing the costs of the training sequences in terms of SM and discretizing to the minimum cost ... ');end                            
                            KM = getUpdatedCosts(Xtrain{sample},params.phmm.SM{k});
                            [~,params.phmm.Dtrain{k}{l}{sample}] = min(KM);
                        elseif ~strcmp(params.phmm.clustType,'none')
                            params.phmm.Dtrain{k}{l}{sample} = Xtrain{sample}';
                        else
                            error('testHMM:noFcn','This discrete case is not considered/implemented');
                        end
                    end
    %                 params.phmm.Dtrain{k}{l} = cell2mat(params.phmm.Dtrain{k}{l});
                else
                    if ~isempty(params.phmm.SM{k})
                         % Obtain Subgesture cost representations for training data
                        display('Computing the costs of the training sequences in terms of SM and discretizing to the minimum cost ... ');
                        KM = getUpdatedCosts(Xtrain,params.phmm.SM{k});
                        [~,params.phmm.Dtrain{k}{l}] = min(KM);
                    elseif ~strcmp(params.phmm.clustType,'none')
                         params.phmm.Dtrain{k}{l} = Xtrain';
                    else
                        error('testHMM:noFcn','This discrete case is not considered/implemented');
                    end
                end
            else
                warning('testHMM:noFcn','Continuous version is not implemented yet');
            end
            %% Learn HMMs
            display(sprintf('Learning the HMM Model for gesture %d ...\n',l));
            if ~params.phmm.pmtk                
                [params.phmm.hmmTR_f{k}{l}, params.phmm.hmmE_f{k}{l}] = ...
                    learnHMM(params.phmm.states,params.phmm.Dtrain{k}{l},params.phmm.it);
            else
                if strcmp(params.phmm.varType,'discrete')
                    [params.phmm.model{k}{l},~] = hmmFit(params.phmm.Dtrain{k}{l}', params.phmm.states, params.phmm.varType);
                elseif strcmp(params.phmm.varType,'gauss') || strcmp(params.phmm.varType,'mixgausstied')
                    for samp = 1:length(Xtrain)
                        Xtrain{samp} = Xtrain{samp}';
                    end
                    io=0;
                    flg=false;
                    while (io<10 && ~flg), %%% try several times
                        try
                            [params.phmm.model{k}{l},~] = hmmFit(Xtrain', params.phmm.states, params.phmm.varType, 'nmix', 2);
                            flg=true;
                        catch e
                            params.phmm.model{k}{l}=[]; io=io+1; %lasterr
                        end
                    end
                    if isempty(params.phmm.model{k}{l})
                        error(e.identifier,e.message);
                    end
                end
            end
        end
        
        %% Evaluating sequences for each model learnt
        for l = 1:length(Xtrain_l)-1
            %% Obtain Training data
            % Xval = Xdev{2}(Ydev{2}.Lfr==l,:);    % gesture vector
            if params.sw > 0
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
                    Xval=Xdev{2}(seg,:); Yval.Lfr=Ydev{2}.Lfr(seg);                         % gesture vector
%                     Xval = Xval(Yval.Lfr == l,:);  % baseline 1+2: get current gesture label
                end
            else
                Yval.Lfr = Ydev{2}.Lfr(Ydev{2}.Lfr==l); Xval = Xval_l{l};                   % gesture cells
            end            
            %% Discretize validation data
            if strcmp(params.phmm.varType,'discrete')
                if ~strcmp(params.phmm.clustType,'none')
                    display('Discretizing validation sequence/s in Key Poses ...');
                    Xval = discretizeSequence(params.phmm.C{k},Xval);
                end
                display('Computing the costs of the validation sequences in terms of SM and discretizing to the minimum cost ... ');
                if iscell(Xval)
                    Dval = cell(1,length(Xval));
                    for sample = 1:length(Xval)
                        if ~isempty(params.phmm.SM{k})
                            %% Obtain Subgesture cost representations for validation data
                            KT = getUpdatedCosts(Xval{sample},params.phmm.SM{k});
                            [~,Dval{sample}] = min(KT);
                        elseif ~strcmp(params.phmm.clustType,'none')
                            Dval{sample} = Xval{sample};
                        else
                            error('testHMM:noFcn','This discrete case is not considered/implemented');
                        end
                    end
                else
                    if ~isempty(params.phmm.SM{k})
                        KT = getUpdatedCosts(Xval,params.phmm.SM{k});
                        [~,Dval] = min(KT);
                    elseif ~strcmp(params.phmm.clustType,'none')
                        Dval = Xval;
                    else
                        error('testHMM:noFcn','This discrete case is not considered/implemented');
                    end
                end
            else
                warning('testHMM:noFcn','Continuous version is not implemented yet');
            end
                
            %% Evaluate sequences and obtain learning rates
            if params.phmm.hmm
                sw = params.sw;
                if sw > 0
                    r = inf;
                    while any(r > length(Dval)-sw)
                        r=randperm(round(length(Dval)),1);
                    end
                    seg=r(1):min(r(1)+sw,length(Yval.Lfr));
                    Dval=Dval(seg);
                    Yval.Lfr=Yval.Lfr(seg);
                end
                [model,score] = evalswHMM(params, Dval, Yval, params.phmm.hmmTR_f{k},params.phmm.hmmE_f{k},params.phmm.model{k});
                return;
            end
            params.phmm.pTrain_f{k}{l} = zeros(length(params.phmm.Dtrain{k}),length(params.phmm.Dtrain{k}{l}));
            params.phmm.mapHMMtrain{k}{l} = zeros(1,length(params.phmm.Dtrain{k}{l}));
            if iscell(Dval), 
                params.phmm.pVal_f{k}{l} = zeros(length(Xtrain_l)-1,length(Dval));
                params.phmm.mapHMMval{k}{l} = zeros(1,length(Dval));
            else
                params.phmm.pVal_f{k}{l} = zeros(length(Xtrain_l)-1,1);                
            end
            
            for m = 1:length(Xtrain_l)-1
                display(sprintf('Evaluation of learning sequences of gesture %d with the HMM model of gesture %d ...',l,m));
                if strcmp(params.phmm.varType,'discrete')
                    if iscell(params.phmm.Dtrain{k}{l})
%                         params.phmm.pTrain_f{k}{l}(m,:) = evaluateSequences([],params.phmm.Dtrain{k}{l},params.phmm.hmmTR_f{k}{m},params.phmm.hmmE_f{k}{m});
                        params.phmm.pTrain_f{k}{l}(m,:) = evalswHMM(params,params.phmm.Dtrain{k}{l},params.phmm.hmmTR_f{k}{m},params.phmm.hmmE_f{k}{m},params.phmm.model{k}{m});
                    end
%                     params.phmm.pVal_f{k}{l}(m,:) = evaluateSequences([],Dval,params.phmm.hmmTR_f{k}{m},params.phmm.hmmE_f{k}{m});
                    params.phmm.pVal_f{k}{l}(m,:) = evalswHMM(params,Dval,params.phmm.hmmTR_f{k}{m},params.phmm.hmmE_f{k}{m},params.phmm.model{k}{m});
                elseif params.phmm.pmtk && strcmp(params.phmm.varType,'gauss') || strcmp(params.phmm.varType,'mixgausstied')
                    params.phmm.pTrain_f{k}{l}(m,:) = evalswHMM(Xtrain_l{l}(1:1+(nSampGest-1)),[],[],params.phmm.model{k}{m});
                    params.phmm.pVal_f{k}{l}(m,:) = evalswHMM(Xval,[],[],params.phmm.model{k}{m});                    
                else
                    error('testHMM:hmmTrainError','Error on the HMM training settings. Check varType and clustType parameters');
                end
            end
            [~,params.phmm.mapHMMtrain{k}{l}] = max(params.phmm.pTrain_f{k}{l});
            if size(params.phmm.pVal_f{k}{l},2) > 1
                [~,params.phmm.mapHMMval{k}{l}] = max(params.phmm.pVal_f{k}{l});
            else
                [~,params.phmm.mapHMMval{k}{l}] = max(params.phmm.pVal_f{k}{l});
            end
            params.phmm.accTrain{k}(l) = sum(params.phmm.mapHMMtrain{k}{l}==l)/length(params.phmm.mapHMMtrain{k}{l});
            params.phmm.accLearn{k}(l) = sum(params.phmm.mapHMMval{k}{l}==l)/length(params.phmm.mapHMMval{k}{l});
            
%             Plot results of the model showing learning and predictive capabilities 
%             plotResults(params.phmm.pTrain_f{k}{mapHMM},params.phmm.pVal_f{k}{mapHMM}, ...
%                params.phmm.hmmE_f{k}{mapHMM},params.phmm.states,k);

            %% minimum probability to define the threshold
%             params.phmm.minProb{k}(l) = min(params.phmm.pVal_f{k}{l}(l,:));
%             hits = sum(params.phmm.pTrain_f{k}{l}(l,:) > params.phmm.minProb{k}(l));
%             params.phmm.accTrain{k}(l) = mean(hits/length(params.phmm.pTrain_f{k}{l}));
%             hits = sum(params.phmm.pVal_f{k}{l}(l,:) > params.phmm.minProb{k}(l));
%             params.phmm.accLearn{k}(l) = mean(hits/length(params.phmm.pVal_f{k}{l})); 
        end
    end
    display(sprintf('Saving model and results ...'));
    save(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'),'params');
    display('Done!');
else    
    load(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'));
    
    %% Choose best fold
    display('Showing Learning results for each fold and gesture ...');
    maxAcc = -inf;
    for k = 1:params.phmm.folds        
        if maxAcc < mean(params.phmm.accLearn{k});
            maxAcc = mean(params.phmm.accLearn{k});
            kBest = k;
        end
    end
    %% Test results
    params.phmm.pTest_f = zeros(1,length(Xtrain_l)-1); params.phmm.accTest = zeros(1,length(Xtrain_l)-1);
    for l = 1:length(Xtrain_l)-1        
%         plotResults(params.phmm.pTrain_f{kBest}{l},params.phmm.pVal_f{kBest}{l},...
%             params.phmm.hmmE_f{kBest}{l},params.phmm.states,kBest);
        if ~exist(strcat('results/',DATATYPE,'/validation/hmm/resProbTestSeqs.mat'),'file'),   
            %% Test                     
            display('Computing the costs of the test sequences in terms of SM and discretizing to the minimum cost ... ');
            Dtest = cell(1,length(Xtest));
            if iscell(Xtest)
                for sample = 1:length(Xtest)
                    KT = getUpdatedCosts(Xtest{sample},params.phmm.SM{kBest});
                    [~,Dtest{sample}] = min(KT);
    %                         Dval{sample} = Xtest{sample};
                end
            else
                KT = getUpdatedCosts(Xtest,params.phmm.SM{kBest});
                [~,Dtest] = min(KT);
    %                     Dtest = Xtest;
            end
            display('Evaluating the final Models with the test set...');
            params.phmm.pTest_f(l) = evaluateSequences([],Dtest,params.phmm.hmmTR_f{kBest}{l},params.phmm.hmmE_f{kBest}{l});
            hits = sum(params.phmm.pTest_f(l) > params.phmm.minProb{kBest}(l));
            params.phmm.accTest(l) = mean(hits/length(params.phmm.pTest_f(l)));

            plotResults(-1,params.phmm.pTest_f(l),params.phmm.hmmE_f{kBest}{l},params.phmm.states,kBest);
        end
    end
    if ~exist(strcat('results/',DATATYPE,'/validation/hmm/resProbTestSeqs.mat'),'file'),   
        save(strcat('results/',DATATYPE,'/validation/hmm/resProbTestSeqs.mat'),'Xtest','params');
    else        
        figure,
        hold on
        plot(params.phmm.minProb{kBest},'black');
        title('Blue: validation set of the selected fold. Red: Test set');
        plot(params.phmm.pTest_f,'blue');
        ylabel('Probability values');
        xlabel('Sequence number');
        hold off
    end
end
