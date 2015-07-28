function s = fitnessFcn(I,X,XtrainT,Xtrain_l,Ytrain,Xval_l,Xval,Yval,params)
    % Validate the data sequences by means of either k-means-DTW and mean DTW models.

    % output:
    %   s: score
    % input:
    %   I: Individual for the genetic algorithm
    %   X: whole training data    
    %   XtrainT: the training data sequence for testing (without noise)
    %   Xtrain_l: training data grouped (labeled) by classes (gestures)
    %   Ytrain: training data labels    
    %   Xval: the validation data sequence for testing (including noise)    
    %   Xval_l: validation data grouped by classes (gestures)
    %   Yval: validation data labels
    %   params: struct with the required parameters to construct the model
    %      N: Number of segments to split the training sequence
    %      N0: Number of segments to split the model sequence
    %      nmin: minimum subsequence width
    %      nmax: maximum sequence width
    %      versions: string with the versions of the kmeans DTW algorithm
    %      dist: distance metric for the k-means DTW
    %      k0: initial data clusters 
    %      M: mean sequence models
    %      bestThs: best thresholds for the k detected gestures    
    %      nThreshs: Number of thresholds for testing
    %      params.D: Dissimilarity matrix
    %      phmm: struct with the required parameters to learn hmm    
    %   state: GA state of the current generation

    global CACHE; global PREDICTIONS; global BESTIND;
    
    [I2,idxUniques] = unique(I,'rows');    
    idxRepsI = ~ismember(1:size(I,1),idxUniques);
    [~,posRepsI2] = ismember(I(idxRepsI,:),I2,'rows');
    if ~isempty(idxRepsI)
        if ~all(isequal(I2(posRepsI2,:),I(idxRepsI,:)))
            error('fitnessFcn:uniqueFailed','All repeated individuals must be members of the list of unique individuals.\nCheck the use of unique');
        end
    end   
    [valid,err] = validateI(I2,params,size(X,1),X);
    % create here new random individuals for the non-valid ones & restore
    % them err = 0. Use the idea of function createInitPopul, and
    % generalize to the case of vectorized on|off (full population in I | 
    % just one individual in I)
    
    [exists,s2] = getCacheVal(int32(round(I2)),params);
    predictions = cell(1,size(I2,1));
    
    if strcmp(params.scoreMeasure,'levenshtein') && any(err)
        err(err < 0) = -err(err < 0) + 21;
    end
    s(idxUniques) = err;
    s(idxRepsI) = err(posRepsI2);
    
    if ~any(valid) || all(exists)
        if all(exists)
            s(idxUniques)=s2;
            s(idxRepsI) = s2(posRepsI2);
        end
        return;
    end
    clear s;
    
    idx = exists==0 & valid==1;
    [I2,k,seg,mk,mnsegs] = decode(I2(idx,:),params);
%     cd('results/temp/');
%     save('temp.mat','seg','X','XtrainT','Ytrain','params','k');
    if params.classification, Xv = Xval_l; else Xv = Xval; end
    if length(idx) > 1
        sc = zeros(1,length(k));
        if params.classification, sc2 = zeros(length(k),3); else sc2 = zeros(length(k),4); end
        preds = cell(1,length(k));
        model = cell(1,length(k));
        % 1) Launch length(k) processes.
        %   Each process calls to the evaluation function for that k,seg,
        %   and it saves a file in a directory with the trained model, 
        %   scores, and predictions
        % 2) Save the results of each process in the structures model{i}, sc(i), preds{i}
        % 3) Release memory, removing the files of the evaluation
%         f = zeros(1,length(k));
%         for i = 1:length(k)
%             cmd = sprintf('#!/bin/bash\nqsub -t %i -q short.q -l mem=2G subgesture_eval.sh',i);
%             f(i) = system(cmd);
%         end
%         for i = 1:length(k)
%             waitfor(f(i));
%             if exist(strcat('eval',num2str(i),'.mat'),'file')
%                 load(strcat('eval',num2str(i),'.mat'));
%                 delete(strcat('eval',num2str(i),'.mat'));
%             else
%                 error('fitnessFcn:FileNotFound','File for the evaluation not found');
%             end            
%             model{i} = m; sc(i) = sco; preds{i} = predic; 
%             clear m sco predic;            
%         end
%         if exist('temp.mat','file'), delete('temp.mat'); end;
        parfor i = 1:length(k)            
            if iscell(seg)
                segA = seg{i};
            else
                segA = reshape(seg(i,:,:),size(seg,2),size(seg,3));
            end
            model{i} = evalFit(X,XtrainT,Xtrain_l,Ytrain,params,k(i),segA,mnsegs(i),mk(i));
            Yv = Yval;      % Parfor doesn't accept modifying the value of Yval directly
            display('Optimizing model parameters over validation ...');
            if ~params.phmm.hmm
%                 model{i}.sw = 0;           % Evaluate the whole validation sequence
                if params.darwin
                    display('Obtaining Validation Sequences from Subgesture Ranks ...')
                    [~,model{i}.KT] = computeSGFromU(model{i}.Us,Xv);
                end
                if params.svm
                    sc(i) = testDarwin(model{i}.KM, model{i}.KT);
                else
                    [model{i},sc(i),~,preds{i}] = g(model{i},Xv,Yval);    % learn&optimize over validation
                end
            else
                if ~strcmp(params.phmm.clustType,'none')
                    display('Discretizing validation sequences in Key Poses ...');
                    Xva = discretizeSequence(model{i}.Ctrain,Xv);      
                else
                    Xva = Xv;
                end
                %% Discretize validation sequence
                display('Discretizing validation sequences in Subgestures ...');
                sw = model{i}.sw;
                if sw > 0 && ~params.classification
                    r = inf;
                    while any(r > length(Xva)-sw)
                        r=randperm(round(length(Xva)),1);
                    end
                    sseg=r(1):min(r(1)+sw,length(Yv.Lfr));
                    Xva=Xva(sseg,:);
                    Yv.Lfr=Yv.Lfr(sseg);
                end
                if iscell(Xva)
                    if params.darwin
                        Dval = computeSGFromU(model{i}.Us,Xva);
                    else
                        Dval = cell(1,length(Xva));
                        for l = 1:length(Xva)
                            Dval{l} = cell(1,length(Xva{l}));
                            for sample = 1:length(Xva{l})
                                KT = getUpdatedCosts(Xva{l}{sample},model{i}.SM);
                                [~,Dval{l}{sample}] = min(KT);
                            end
                        end
                    end
                else
                    if params.darwin
                        Dval = computeSGFromU(model{i}.Us,Xva);
                    else
                        KT = getUpdatedCosts(Xva,model{i}.SM);
                        [~,Dval] = min(KT);
                    end
                end
%                 sc(i) = evaluateHMM(Dval, model{i}.phmm.hmmTR, model{i}.phmm.hmmE);
                [model{i},sc(i),sc2(i,:)] = evalswHMM(model{i}, Dval, Yval)
                preds{i} = 0;
            end
        end
        s2(idx) = sc;
        predictions(idx) = preds;
        clear preds;
    else
        predictions = cell(1);
        model = cell(1);
        if iscell(seg)
            model{1} = evalFit(X,XtrainT,Xtrain_l,Ytrain,params,k,cell2mat(seg),mnsegs,mk);
        else
            model{1} = evalFit(X,XtrainT,Xtrain_l,Ytrain,params,k,reshape(seg,size(seg,2),size(seg,3)),mnsegs,mk);
        end
        display('Optimizing model parameters over validation ...');
        if ~params.phmm.hmm
%             model{1}.sw = 0;           % Evaluate the whole validation sequence 
            if params.darwin
                display('Obtaining Validation Sequences from Subgesture Ranks U ...')
                [~,model{1}.KT] = computeSGFromU(model{1}.Us,Xv);
            end
            if params.svm
                s2 = testDarwin(model{1}.KM, model{1}.KT); sc2 = -inf;
            else
                [model{1},s2,sc2,predictions{1}] = g(model{1},Xv,Yval);   % learn&optimize over validation
            end
            s2
        else
            if ~strcmp(params.phmm.clustType,'none')
                display('Discretizing validation sequences in Key Poses ...');
                Xva = discretizeSequence(model{1}.Ctrain,Xv);
            else
                Xva = Xv;
            end
            %% Discretize validation sequence
            display('Discretizing validation sequences in Subgestures ...');
            sw = model{1}.sw;
            if sw > 0 && ~params.classification
                r = inf;
                while any(r > length(Xva)-sw)
                    r=randperm(round(length(Xva)),1);
                end
                sseg=r(1):min(r(1)+sw,length(Yval.Lfr));
                Xva=Xva(sseg,:);
                Yval.Lfr=Yval.Lfr(sseg);
            end
            if iscell(Xva)
                if params.darwin
                    Dval = computeSGFromU(model{1}.Us,Xva);
                else
                    Dval = cell(1,length(Xva));
                    for l = 1:length(Xva)
                        Dval{l} = cell(1,length(Xva{l}));
                        for sample = 1:length(Xva{l})
                            KT = getUpdatedCosts(Xva{l}{sample},model{1}.SM);
                            [~,Dval{l}{sample}] = min(KT);
                        end
                    end
                end
            else
                if params.darwin
                    Dval = computeSGFromU(model{1}.Us,Xva);
                else
                    KT = getUpdatedCosts(Xva,model{1}.SM);
                    [~,Dval] = min(KT);
                end
            end            
            [model{1},s2,sc2] = evalswHMM(model{1}, Dval, Yval);
            predictions{1} = 0;
            s2
        end
    end

    %%%% if want to plot the performance of the best solution so far
%     plotmistakes(predictions{best}, Yval,0)

     % Add to cache (FILO stack mode avoiding to overwrite the highest value)
    pos = CACHE.pos;
    posEnd = pos+sum(idx)-1;
    remain = 0;
    if posEnd >= pos
        if length(CACHE.eval) < posEnd
            remain = posEnd - length(CACHE.eval);
            posEnd = length(CACHE.eval);
        end
        idxEnd = posEnd-pos+1;
        idxC = pos:posEnd;
        [maxSC,i] = max(CACHE.eval);
        if maxSC > 0
            [~,loc] = ismember(i,idxC);
            if any(loc)
                idxC(loc:end) = idxC(loc:end) + 1;
                addRemain = idxC > length(CACHE.eval);
                if any(addRemain)
                    idxC(addRemain) = 1:sum(addRemain);
                end
            end
        end        
        CACHE.ind(idxC,:) = I2(1:idxEnd,:); 
        s3 = s2(idx); 
        CACHE.eval(idxC) = s3(1:idxEnd);
        if all(isinf(sc2))
            s4 = sc2(idx,:);
            CACHE.scores(idxC,:) = s4(1:idxEnd,:);
        end
        if remain
            pos = 1;
            posEnd = pos+remain-1;
            idxC = pos:posEnd;
            [maxSC,i] = max(CACHE.eval);
            if maxSC > 0
                [~,loc] = ismember(i,idxC);
                if any(loc)
                    idxC(loc:end) = idxC(loc:end) + 1;
                    addRemain = idxC > length(CACHE.eval);
                    if any(addRemain)
                        idxC(addRemain) = 1:sum(addRemain);
                    end
                end
            end            
            CACHE.ind(idxC,:) = I2(idxEnd+1:end,:); 
            CACHE.eval(idxC) = s3(idxEnd+1:end);
        end
        CACHE.pos = posEnd + 1;
    end
    
    % Save current best model    
    [bestScore,best] = max(s2(idx));
    if bestScore >= max(CACHE.eval)
        BESTIND(end).model = model{best};
    end
    if isempty(BESTIND(end).model)
        BESTIND(end).model = BESTIND(end-1).model;
    end
    % Return all final scores    
    s(idxUniques) = s2;
    s(idxRepsI) = s2(posRepsI2);
    preds(idxUniques) = predictions;
    preds(idxRepsI) = predictions(posRepsI2);    
    [~,bestPos] = max(s);
    PREDICTIONS = preds{bestPos};
end

function [I2,k,seg,mk,mnsegs] = decode(I,params)
    k = round(I(:,1));
    I2 = zeros(size(I));
    I2(:,1) = k(:);     % get k
    if params.probSeg > 0
        seg = cell(1,size(I,1));
    end
    for i = 1:size(I,1)
        Iseg = round(I(i,2:params.N*2+1));
        I2(i,2:params.N*2+1) = Iseg;
        Iseg = Iseg(Iseg < inf);    % choose non-infinity segments
        if params.probSeg > 0
            seg{i} = round(reshape(Iseg,[2 size(Iseg,2)/2]));        
        end
    end
    if params.probSeg == 0
        seg = round(reshape(I(:,2:params.N*2+1),[size(I,1) 2 (params.N*2)/2]));
    end
    
    mnsegs = zeros(1,size(I,1)); mk = zeros(1,size(I,1));
    % get model segments and k
    global MEDIANTYPE;
    if strcmp(params.mType,MEDIANTYPE{2}) || strcmp(params.mType,MEDIANTYPE{4})
        mnsegs = round(I(:,end-1));
        mk = round(I(:,end));
        I2(:,end-1) = mnsegs;
        I2(:,end) = mk;
    end
end

function [valid,err] = validateI(I,params,maxSeg,X)
    global BASELINE;
    [I2,k,seg,mk,mnsegs] = decode(I,params);
    valid = ones(1,length(k)); 
    
    %% check k for subgestures
    nsegs = sum(transpose(I2(:,2:params.N*2+1) < inf));
    nsegs(mod(nsegs,2) > 0) = nsegs(mod(nsegs,2) > 0) - 1;
    nsegs = nsegs/2;
    kBad = k < params.k0 | k > params.N | k(:) > nsegs(:);
    e0 = zeros(1,length(k));
    if ~isempty(k(kBad))
        e0(kBad) = abs(max(k(kBad)-params.k0,params.N-k(kBad)))/10000;
        valid(kBad) = 0;
    end
    
    %% Check subgesture segments
    e1 = zeros(1,size(I2,1)); e2 = zeros(1,size(I2,1));
    e3 = zeros(1,size(I2,1)); e4 = zeros(1,size(I2,1));
    for i = 1:size(I2,1)
        if iscell(seg)
            e1(i) = sum(seg{i}(1,:) < 1);
        else
            e1(i) = sum(seg(i,1,:) < 1);
        end
        if strcmp(params.Baseline,BASELINE{1})
            if iscell(seg)
                e2(i) = sum(seg{i}(1,:) > maxSeg-params.nmax);
                e3(i) = sum(seg{i}(2,:) < params.nmin);
                e4(i) = sum(seg{i}(2,:) > params.nmax);
            else
                e2(i) = sum(seg(i,1,:) > maxSeg-params.nmax);
                e3(i) = sum(seg(i,2,:) < params.nmin);
                e4(i) = sum(seg(i,2,:) > params.nmax);
            end
        else
            e2(i) = 0;
            if iscell(seg)
                e3(i) = sum(seg{i}(2,:) < 2);
                e4(i) = sum(seg{i}(1,:)+seg{i}(2,:)-1 > maxSeg);
            else
                e3(i) = sum(seg(i,2,:) < 2);
                e4(i) = sum(seg(i,1,:)+seg(i,2,:)-1 > maxSeg);
            end
        end            
    end

    % check non-zero segments
    for i = 1:size(I2,1)
        if iscell(seg)
            l = length(seg{i})-1;
        else
            l = params.N-1;
        end
        for j = 1:l
            if ~any([e1(i) e2(i) e3(i) e4(i)])
                if iscell(seg)
                    in = seg{i}(1,j);
                    fi = in + seg{i}(2,j) - 1;
                else
                    in = seg(i,1,j);
                    fi = in + seg(i,2,j) - 1;
                end
                if ~any(any(X(in:fi,:)))
                    e1(i) = size(I2,1);
                end
            end
        end
    end
        
    %% check number of segments and k for the models    
    e5 = zeros(1,size(I2,1)); 
    e6 = zeros(1,size(I2,1));
    global MEDIANTYPE;
    if strcmp(params.mType,MEDIANTYPE{2}) || strcmp(params.mType,MEDIANTYPE{4})
        nsBad = mnsegs < params.k0 | mnsegs > params.N0;
        if ~isempty(nsegs(nsBad))
            e5(nsBad) = abs(max(mnsegs(nsBad)-params.k0,params.N0-mnsegs(nsBad)))/100;
            valid(nsBad) = 0;
        end
        kBad = mk < params.k0-1 | mk > params.N0 | mk(:) > mnsegs(:);
        if ~isempty(mk(kBad))
            e6(kBad) = abs(max(mk(kBad)-params.k0-1,params.N0-mk(kBad)))/100;
            valid(kBad) = 0;
        end
    end
    
    %% calculate final error
    e = e1+e2+e3+e4+e5;
    if any(e)
        valid(e~=0) = 0;
        % If there is more than the half population invalid, create new
        % random initial population.
    end
    err = - (0.25*e0 + 0.5*e/length(I2(:,2:end))/2 + 0.25*e6);
    if err, 
        warning('fitnessFcn:invInd','Invalid individuals detected in Fitness function'); 
    end
end

function segs = getDataPartitions(data,partitions)
    segs = cell(1,size(partitions,1));
    for i = 1:length(segs)
        segs{i} = data(partitions(i,1):partitions(i,1)+partitions(i,2)-1,:);
    end
end

function [exists,s] = getCacheVal(I,params)
    global CACHE;
    
    if strcmp(params.scoreMeasure,'overlap')
        s=-1*ones(1,size(I,1));
    elseif strcmp(params.scoreMeasure,'levenshtein')
        s=inf*ones(1,size(I,1));
    end
    exists = zeros(1,size(I,1));        
    if ~isempty(CACHE.ind)
%         tic;
        [members,loc] = ismember(I,CACHE.ind,'rows');
%         toc;
        if any(members)
            exists(members) = 1;             
            s(members) = CACHE.eval(loc(loc > 0));
            if strcmp(params.scoreMeasure,'overlap') && ...
                    ~ismember(max(CACHE.eval),s) && length(s) > 1 || ...
                    strcmp(params.scoreMeasure,'levenshtein') && ...
                    ~ismember(min(CACHE.eval),s) && length(s) > 1
                warning('fitnessFcn:ElitismError','Elitist members not found within the population.\nCurrent position: %d\nBest in cache: %.4f\nBest in s (elitist): %.4f',CACHE.pos,max(CACHE.eval),max(s));
            end
        end
    end
end

function model = evalFit(X,XtrainT,Xtrain_l,Ytrain,params,k,seg,mnseg,mk)
% Output:
%   s: scores
%   model: model structure
%      SM: Subgesture models (centroid sequences) for the k gestures
%      C: cluster indices for each sample
%      D: Dissimilarity matrix between subgesture models SM
%      M: mean sequence models
%      KM: Cost matrix of the mean sequence models in terms of SM
%      KT: Cost matrix of the test sequence in terms of SM    
%      bestThs: best thresholds for the k detected gestures    
%      nThreshs: Number of thresholds for testing

    % Obtain clsuters using the k-means DTW algorithm    
    if any(k > size(seg,2))
        error('fitnessFcn:k','k cannot be greater than the number of segments');
    end
    
    X_I = getDataPartitions(X,seg');    
    model = params;
    
    if params.darwin
        %% Computing Video Darwin Subgesture Us
        display('Computing Subgesture Ranks U for each X partition ...')
        for i=2:length(X_I)+1
            W = genRepresentation(X_I{i-1},1);
%             if mod(i,10)==0, display(sprintf('Current gesture segment %d of %d\n',i-1,length(X_I))); end
            Xtrain(i-1,:)=W';
        end
        display('Done!');
        [~,model.Us] = kmeans(Xtrain,k,'EmptyAction','singleton');
        model.D = getSimilarities(model.Us,model);
        if sum(~any(model.D))
            error('fitnessFcn:D','Some elements of the dissimilarity matrix are wrong');
        end
    else
        %% Temporal Clustering
        display('Computing Subgestures as Temporal Clusters for each X partition ...')
        [CsTrain,~,mErrsV,~,timeV,~,Z] = runKMeansDTW(params,k,k,[],[],[],Ytrain,[],X_I,[]);
        [~,kV] = min(mErrsV);
        model.C = CsTrain{kV}{timeV(kV)};

        %% Obtain Subgesture Model for training/learning data
        model.SM = Z{kV}{timeV}; emptyCells = cellfun(@isempty,model.SM); model.SM(emptyCells) = [];
    end
    if ~params.phmm.hmm
        if params.darwin
            display('Obtaining Subgesture Models from Subgesture Ranks U ...')
            if params.svm
                [~,model.KM] = computeSGFromU(model.Us,Xtrain_l);
            else
                model.M = cell(1,length(params.M)); model.KM = cell(1,length(params.M)); 
                for i = 1:length(params.M)
                    [~,model.KM{i}] = computeSGFromU(model.Us,params.M{i});
                end
                model.M = params.M;
            end
        else
            %% compute Similarity matrix
            if ~params.pdtw
                model.D = getSimilarities(model.SM,model);
                if sum(~any(model.D))
                    error('fitnessFcn:D','Some elements of the dissimilarity matrix are wrong');
                end
            end

            %% Compute Median/Max (SubGesture) Models for each gesture 
            if (strcmp(params.mType,'modelSM1') || strcmp(params.mType,'allSM1')) ... 
                    && mnseg > 0 && mk > 0
                model.M = getSM(params,Xtrain_l,model,mnseg,mk);
            elseif strcmp(params.mType,'modelSM2') || strcmp(params.mType,'allSM2')
                model.M = getSM(params,Xtrain_l,model);
            elseif strcmp(params.mType,'direct') || strcmp(params.mType,'DCSR')
                model.M = params.M;
            else
                error('fitnessFcn:optError','Option chosen for the Subgesture Models is incorrect. Check the value of params.mType');
            end
            %% Compute costs of representing (Subgesture) Models 'M' in terms of Subgesture Models 'SM'
            model.KM = cell(1,length(Xtrain_l));
            if strcmp(params.mCostType,'mean')
                display('Computing the costs of the training sequences in terms of SM and mean align to the minimum cost ... ');
                for l = 1:length(model.KM)
                    model.KM{l} = cell(1,length(Xtrain_l{l}));
                    for sample = 1:length(model.KM{l})
                        KM = getUpdatedCosts(Xtrain_l{l}{sample},model.SM);
                        model.KM{l}{sample} = KM';
    %                     model.KM{l}{sample} = min(KM)';
                    end
                end
                model.KM = getModels(model.KM,length(model.KM),params);
                for l = 1:length(model.KM)
                    model.KM{l} = model.KM{l}';
                end
            else
    %             tic;
                display('Computing the costs of the models M in terms of SM ...');
                for i = 1:length(model.KM)            
                    if params.k > 0
                        model.KM{i} = getUpdatedCosts(model.M{i}{params.k},model.SM);
                    else
                        model.KM{i} = getUpdatedCosts(model.M{i},model.SM);
                    end
                end
    %             toc;
            end
        end
    else
        model.phmm.hmmTR_f = cell(1,length(Xtrain_l)); 
        model.phmm.hmmE_f = cell(1,length(Xtrain_l));
        model.phmm.model = cell(1,length(Xtrain_l));

        if params.darwin, 
            display('Obtaining Subgesture Models from Subgesture Ranks U ...')
            MD = computeSGFromU(model.Us,Xtrain_l); 
        end
        display('Discretizing training sequence and learning the HMM Models for each gesture');
        for l = 1:length(Xtrain_l)
            if params.darwin
                Dtrain = MD{l};
            else
                Xtrain = Xtrain_l{l};
                if ~strcmp(params.phmm.clustType,'none')
                    %% Discretize gestures key poses
                    model.Ctrain = performClustering(Xtrain,Ytrain,params.phmm.clustType,params.phmm.kD,params.phmm.cIters);
                    Xtrain = discretizeSequence(model.Ctrain,Xtrain);
                end            
                Dtrain = cell(1,length(Xtrain));
                for sample = 1:length(Xtrain)
                    %% Discretize in subgestures
                    KM = getUpdatedCosts(Xtrain{sample},model.SM);
                    [~,Dtrain{sample}] = min(KM);
                end
            end
            %% train HMM for this gesture
            if ~params.phmm.pmtk
                [model.phmm.hmmTR_f{l}, model.phmm.hmmE_f{l}] = learnHMM(params.phmm.states,Dtrain,params.phmm.it);
                model.phmm.model{l} = [];
            else
                [model.phmm.model{l},~] = hmmFit(Dtrain', params.phmm.states, params.phmm.varType);
            end
        end
    end
end