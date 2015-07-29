function [score,model,bestScores] = testHMM(params)

% clear all
% close all
% addPath

%% Generate global variables, cache and parameters
varload
if ~exist('params','var')
    params.scoreMeasure = 'overlap';  % Score Measure: 'overlap' or 'levenshtein'    
    clear CACHE S BASELINE OPTIONS STATE;
end
if ~isfield('params','scoreMeasure'),
    params.scoreMeasure = 'overlap';  % Score Measure: 'overlap' or 'levenshtein'    
end
nSampGest = 10;
switch params.score2optim
    case 'o', if ~params.classification, params.score2optim = 1; else error('g:optErr','Spotting is not allowed in classification ...for now... (ODL!)'); end
    case 'p', if ~params.classification, params.score2optim = 2; else params.score2optim = 1; end
    case 'r', if ~params.classification, params.score2optim = 3; else params.score2optim = 2; end
    case 'a', if ~params.classification, params.score2optim = 4; else params.score2optim = 3; end
end

%% Prepare training data depending on the chosen option and parameters
% Load data:
%     if nframesSeg is 0, then initial segmentation is generated from the skeleton labels
COORDS = 'world'; NAT = 3;
[X,Y,Xtest,Ytest] = prepareData(nrsamples,nseqs,nframesSeg,params.k0);
% display('Press a key to continue...');
% pause();

%% Obtain all training samples grouped (labeled) by gestures
Xtrain_l = getGroupedGestures(X,Y,1); if sum(cellfun(@isempty,Xtrain_l)), error('Empty gesture classes'); end
Xval_l = getGroupedGestures(X,Y,2); if sum(cellfun(@isempty,Xval_l)), error('Empty gesture classes'); end

%% Generate learning sequences
% l = [24 78 150];    % 78 (more samples for each gesture when k=3);
l = [];
[Xdev,Ydev] = getDevSequences(X,Y,l,noise,secsBatch,0);
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
% params.phmm.pTrain_f = cell(1,params.phmm.folds); params.phmm.pVal_f = cell(1,params.phmm.folds); 
params.phmm.predTrain = cell(1,params.phmm.folds); params.phmm.predVal = cell(1,params.phmm.folds); 
% params.phmm.mapHMMtrain = cell(1,params.phmm.folds); params.phmm.mapHMMval = cell(1,params.phmm.folds); params.phmm.minProb = cell(1,params.phmm.folds); 
% params.phmm.accTrain = cell(1,params.phmm.folds); params.phmm.accLearn = cell(1,params.phmm.folds);

%% Train and save learning results
if ~exist(strcat('results/',DATATYPE,'/validation/hmm/learningResults.mat'),'file'),
    
    for k = 1:params.phmm.folds,
        display(sprintf('\n Fold %d',k));
        if strcmp(params.phmm.varType,'discrete') 
            params.phmm.Dtrain{k} = cell(1,length(Xtrain_l)-1);
        end
%         params.phmm.pTrain_f{k} = cell(1,length(Xtrain_l)-1); params.phmm.pVal_f{k} = cell(1,length(Xval_l)-1);        
%         params.phmm.mapHMMtrain{k} = cell(1,length(Xtrain_l)-1); params.phmm.mapHMMval{k} = cell(1,length(Xval_l)-1);
%         params.phmm.minProb{k} = zeros(1,length(Xtrain_l)-1);
%         params.phmm.accTrain{k} = cell(1,length(Xtrain_l)-1); params.phmm.accLearn{k} = cell(1,length(Xval_l)-1);
        params.phmm.predTrain{k} = cell(1,length(Xtrain_l)-1); params.phmm.predVal{k} = cell(1,length(Xval_l)-1);
        params.phmm.hmmTR_f{k} = cell(1,length(Xtrain_l)-1); params.phmm.hmmE_f{k} = cell(1,length(Xval_l)-1);
        params.phmm.model{k} = cell(1,length(Xtrain_l)-1); params.phmm.path{k} = cell(1,length(Xval_l)-1);
        
        %% Select a sample containing gesture samples of each class
        Xtrain = [];
        for l = 1:length(Xtrain_l)
            if ~nSampGest, 
                Xtrain = [Xtrain Xtrain_l{l}];
            else
                Xtrain = [Xtrain Xtrain_l{l}(1:min(length(Xtrain_l{l}),nSampGest))];
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
        Ytrain = cell(1,length(Xtrain_l)-1);
        for l = 1:length(Xtrain_l)-1
            %% Obtain Training data for each gesture class
            if ~nSampGest
                Xtrain = Xtrain_l{l};
            else
                Xtrain = Xtrain_l{l}(1:min(length(Xtrain_l{l}),nSampGest));
%                 Xtrain = cell(1,length(Xtrain_l{l}));
%                 for s = 1:length(Xtrain_l{l})
%                     Xtrain{s} = Xtrain_l{l}{s}(3:end-2,:);
%                 end
            end
            len = cellfun(@length,Xtrain); Ytrain{l} = cell(1,length(Xtrain));
            for g = 1:length(len), Ytrain{l}{g} = l*ones(1,len(g)); end
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
                    [params.phmm.model{k}{l},~] = hmmFit(params.phmm.Dtrain{k}{l}', params.phmm.states, params.phmm.varType,'verbose',false);
                elseif strcmp(params.phmm.varType,'gauss') || strcmp(params.phmm.varType,'mixgausstied')
                    for samp = 1:length(Xtrain)
                        Xtrain{samp} = Xtrain{samp}';
                    end
                    io=0;
                    flg=false;
                    while (io<10 && ~flg), %%% try several times
                        try
                            [params.phmm.model{k}{l},~] = hmmFit(Xtrain', params.phmm.states, params.phmm.varType, 'nmix', 2,'verbose',false);
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
                if params.phmm.hmm
                    Xval = Xdev{2}; Yval.Lfr = Ydev{2}.Lfr;
                else
                    if ~nSampGest
                        Xval = Xval_l{l};
                    else
                        Xval = Xval_l{l}(1:min(length(Xval_l{l}),nSampGest));
%                         Xval = cell(1,length(Xval_l{l}));
%                         for s = 1:length(Xval_l{l})
%                             Xval{s} = Xval_l{l}{s}(3:end-2,:);
%                         end
                    end
                    len = cellfun(@length,Xval); Yval = cell(1,length(Xval));
                    for g = 1:length(len), Yval{g} = l*ones(1,len(g)); end
                end
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
                [model,score,bestScores] = evalswHMM(params, Dval, Yval);
                model.SM = params.phmm.SM{k};
                return;
            end
            params.phmm.predTrain{k}{l} = zeros(1,length(params.phmm.Dtrain{k}{l}));
%             params.phmm.mapHMMtrain{k}{l} = zeros(1,length(params.phmm.Dtrain{k}{l}));
            if iscell(Dval), 
                params.phmm.predVal{k}{l} = zeros(1,length(Dval));
%                 params.phmm.mapHMMval{k}{l} = zeros(1,length(Dval));
            end
            
            display(sprintf('Evaluation of learning sequences of gesture %d with the HMM models ...',l));
            if strcmp(params.phmm.varType,'discrete')
                if iscell(params.phmm.Dtrain{k}{l})
%                         params.phmm.pTrain_f{k}{l}(m,:) = evaluateSequences([],params.phmm.Dtrain{k}{l},params.phmm.hmmTR_f{k}{m},params.phmm.hmmE_f{k}{m});
                    [~,params.phmm.predTrain{k}{l},~] = evalswHMM(params,params.phmm.Dtrain{k}{l},Ytrain{l});
                end
%                     params.phmm.pVal_f{k}{l}(m,:) = evaluateSequences([],Dval,params.phmm.hmmTR_f{k}{m},params.phmm.hmmE_f{k}{m});
                [~,params.phmm.predVal{k}{l},~] = evalswHMM(params, Dval, Yval);
            elseif params.phmm.pmtk && strcmp(params.phmm.varType,'gauss') || strcmp(params.phmm.varType,'mixgausstied')
                [~,params.phmm.predTrain{k}{l},~] = evalswHMM(params, Xtrain_l{l}(1:1+(nSampGest-1)), Ytrain);
                [~,params.phmm.predVal{k}{l},~] = evalswHMM(params, Xval, Yval);
            else
                error('testHMM:hmmTrainError','Error on the HMM training settings. Check varType and clustType parameters');
            end            
%             [~,params.phmm.mapHMMtrain{k}{l}] = max(params.phmm.pTrain_f{k}{l});
%             if size(params.phmm.pVal_f{k}{l},2) > 1
%                 [~,params.phmm.mapHMMval{k}{l}(m)] = max(params.phmm.accTest{k}{l});
%             else
%                 [~,params.phmm.mapHMMval{k}{l}] = max(params.phmm.accTest{k}{l});
%             end

%             Plot results of the model showing learning and predictive capabilities 
%             plotResults(params.phmm.pTrain_f{k}{mapHMM},params.phmm.pVal_f{k}{mapHMM}, ...
%                params.phmm.hmmE_f{k}{mapHMM},params.phmm.states,k);
        end
        CFTrain = confusionmat(cell2mat(params.phmm.predTrain{k}),1:length(Xtrain_l)-1);
        imshow(CFTrain)
        CFTest = confusionmat(cell2mat(params.phmm.predVal{k}),1:length(Xval_l)-1);
        imshow(CFTest)
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
