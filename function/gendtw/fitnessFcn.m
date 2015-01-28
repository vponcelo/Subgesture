function s = fitnessFcn(I,Xtrain,XtrainT,Ytrain,Xval,Yval,params)
    % Validate the data sequences by means of either k-means-DTW and mean DTW models.

    % output:
    %   s: score
    % input:
    %   I: Individual for the genetic algorithm
    %   Xtrain: training data
    %   XtrainT: a training data sequence for testing
    %   Xval: validation data
    %   Ytrain: training data labels
    %   Yval: validation data labels
    %   params: struct with the required parameters to construct the model
    %      N: Number of segments to split the training sequence
    %      nmin: minimum subsequence width
    %      nmax: maximum sequence width
    %      versions: string with the versions of the kmeans DTW algorithm
    %      dist: distance metric for the k-means DTW
    %      k0: initial data clusters 
    %      M: mean sequence models
    %      bestThs: best thresholds for the k detected gestures    
    %      nThreshs: Number of thresholds for testing
    %      params.D: Dissimilarity matrix
    %   state: GA state of the current generation

    global CACHE;
    global PREDICTIONS;
    global MODEL;
    
    [I2,idxUniques] = unique(I,'rows');    
    idxRepsI = ~ismember(1:size(I,1),idxUniques);
    [~,posRepsI2] = ismember(I(idxRepsI,:),I2,'rows');
    if ~isempty(idxRepsI)
        if ~all(isequal(I2(posRepsI2,:),I(idxRepsI,:)))
            error('fitnessFcn:uniqueFailed','All repeated individuals must be members of the list of unique individuals.\nCheck the use of unique');
        end
    end   
    [valid,err] = validateI(I2,params,size(Xtrain,1),Xtrain);
    [exists,s2] = getCacheVal(int32(round(I2)),params);
    predictions = cell(1,size(I2,1));
    
    if strcmp(params.scoreMeasure,'levenshtein') && any(err)
        err(err < 0) = -err(err < 0) + 21;
    end
    s(idxUniques) = err;
    s(idxRepsI) = err(posRepsI2);
    
    if ~any(valid) || all(exists)
        return;   
    end
    clear s;
    
    idx = exists==0 & valid==1;
    [I2,k,seg] = decode(I2(idx,:),params);
    if length(idx) > 1
        sc = zeros(1,length(k));
        preds = cell(1,length(k));
        model = cell(1,length(k));
%         tic;
        parfor i = 1:length(k)            
            if iscell(seg)
                segA = seg{i};
            else
                segA = reshape(seg(i,:,:),size(seg,2),size(seg,3));                
            end
            [~,model{i}] = evalFit(Xtrain,XtrainT,Ytrain,params,k(i),segA);
            display('Validating the model ...');
            [~,sc(i),preds{i}] = g(model{i},Xval,Yval);
            display('Done!');
        end
%         toc;
        s2(idx) = sc;
        predictions(idx) = preds;
        clear preds;
    else
        predictions = cell(1);
        model = cell(1);
%         tic;
        if iscell(seg)
            [~,model{1}] = evalFit(Xtrain,XtrainT,Ytrain,params,k,cell2mat(seg));            
        else
            [~,model{1}] = evalFit(Xtrain,XtrainT,Ytrain,params,k,reshape(seg,size(seg,2),size(seg,3)));
        end
        [~,s2,predictions{1}] = g(model{1},Xval,Yval);
%         toc;
    end

    % Save current best model
    [bestScore,best] = max(s2(idx));
    if bestScore > max(CACHE.eval)
        MODEL = model{best};
    end
    
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
        
    % Return all final scores    
    s(idxUniques) = s2;
    s(idxRepsI) = s2(posRepsI2);
    preds(idxUniques) = predictions;
    preds(idxRepsI) = predictions(posRepsI2);    
    [~,bestPos] = max(s);
    PREDICTIONS = preds{bestPos};
end

function [I2,k,seg] = decode(I,params)
    k = round(I(:,1));
    I2 = zeros(size(I));
    I2(:,1) = k(:);
    if params.probSeg > 0
        seg = cell(1,size(I,1));
    end
    for i = 1:size(I,1)
        Iseg = round(I(i,2:end));
        I2(i,2:end) = Iseg;
        Iseg = Iseg(Iseg < inf);
        if params.probSeg > 0
            seg{i} = round(reshape(Iseg,[2 size(Iseg,2)/2]));        
        end
    end
    if params.probSeg == 0
        seg = round(reshape(I(:,2:end),[size(I,1) 2 (size(I,2)-1)/2]));
    end
end

function [valid,err] = validateI(I,params,maxSeg,X)
    global BASELINE;
    [I2,k,seg] = decode(I,params);
    valid = ones(1,length(k)); 
    
    nsegs = sum(transpose(I2 < inf));
    nsegs(mod(nsegs,2) > 0) = nsegs(mod(nsegs,2) > 0) - 1;
    nsegs = nsegs/2;
    kBad = k < params.k0 | k > params.N | k(:) >= nsegs(:);
    e0 = zeros(1,length(k));
    if ~isempty(k(kBad))
        e0(kBad) = abs(max(k(kBad)-params.k0,params.N-k(kBad)))/100;
        valid(kBad) = 0;
    end
    
    % Check segments
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
                e4(i) = sum(seg{i}(2,:) > params.nmax+1);
            else
                e2(i) = sum(seg(i,1,:) > maxSeg-params.nmax);
                e3(i) = sum(seg(i,2,:) < params.nmin);
                e4(i) = sum(seg(i,2,:) > params.nmax+1);
            end
        else
            e2(i) = 0;
            if iscell(seg)
                e3(i) = sum(seg{i}(2,:) < 2);
                e4(i) = sum(seg{i}(1,:)+seg{i}(2,:)-1 > maxSeg);
            else
                e3(i) = sum(seg(i,2,:) < 2);
                e4(i) = sum(seg(i,1,:)+seg{i}(2,:)-1 > maxSeg);
            end
        end            
    end
    
%     if size(I,1) > 1        
%         for i = 1:size(I,1)%             
%             if strcmp(params.Baseline,BASELINE{1})
%                 e2(i) = sum(seg(i,1,:) > maxSeg-params.nmax);
%                 e3(i) = sum(seg(i,2,:) < params.nmin);
%                 e4(i) = sum(seg(i,2,:) > params.nmax+1);
%             else
%                 e2(i) = 0;
%                 e3(i) = sum(seg(i,2,:) < 2);
%                 e4(i) = sum(seg(i,1,:)+seg{i}(2,:)-1 > maxSeg);
%             end            
%         end
%     else
%         e1 = sum(seg(:,1,:) < 1);
%         if strcmp(params.Baseline,BASELINE{1})
%             e2 = sum(seg(:,1,:) > maxSeg-params.nmax & seg(:,1,:) < inf);
%             e3 = sum(seg(:,2,:) < params.nmin);
%             e4 = sum(seg(:,2,:) > params.nmax+1 & seg(:,2,:) < inf);
%         else
%             e2 = 0;
%             e3 = sum(seg(:,2,:) < 2);
%             e4 = sum(seg(:,1,:)+seg(:,2,:)-1 > maxSeg  & seg(:,2,:) < inf);
%         end        
%     end

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
    
    e = e1+e2+e3+e4;
    if any(e)
        valid(e~=0) = 0;
    end    
    err = - (0.5*e0 + 0.5*e/length(I2(:,2:end))/2);
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
                error('fitnessFcn:CacheError','Cache error. Probably caused because elitist members were not in cache\nCurrent position: %d\nBest in cache: %.4f\nBest in s (elitist): %.4f',CACHE.pos,max(CACHE.eval),max(s));
            end
        end
    end
end

function [s,model] = evalFit(Xtrain,XtrainT,Ytrain,params,k,seg)
% Output:
%   model: model structure
%      SM: Subgesture models (centroid sequences) for the k gestures
%      C: cluster indices for each sample
%      D: Dissimilarity matrix between subgesture models SM
%      M: mean sequence models
%      KM: Cost matrix of the mean sequence models in terms of SM
%      KT: Cost matrix of the test sequence in terms of SM    
%      bestThs: best thresholds for the k detected gestures    
%      nThreshs: Number of thresholds for testing

    X_I = getDataPartitions(Xtrain,seg');
    
    if isempty(params.D)
        %% Temporal Clustering
        % Obtain subsets using the k-means DTW algorithm    
        v = char(params.versions);
        if any(k > size(seg,2))
            error('fitnessFcn:k','k cannot be greater than the number of segments');
        end
        [CsTrain,~,mErrsV,~,timeV,~,Z] = runKMeansDTW(v,k,'dtwCost',k,[],[],[],Ytrain,[],X_I,[]);    

        %% Get clustered training/learning data structures     
        [~,kV] = min(mErrsV);
        model.C = CsTrain{kV}{timeV(kV)};

        %% Obtain Subgesture Model for training/learning data        
        model.SM = Z{kV}{timeV};
        emptyCells = cellfun(@isempty,model.SM);
        model.SM(emptyCells) = [];

        %% compute Similarity matrix
        model.D = getSimilarities(model.SM);
        if sum(~any(model.D))
            error('Some elements of the dissimilarity matrix are wrong');
        end
    else
        model.C = params.C;
        model.SM = params.SM;
        model.D = params.D;
        model.KM = params.KM;
    end
            
    %% Compte costs for gesture Models M with respect to the Subgesture Models model.SM
    model.M = params.M;
    model.KM = cell(1,length(params.M));
    display('Computing the costs of the gesture models in terms of SM ...');
    tic;
    for i = 1:length(model.KM)
        if params.k > 0
            model.KM{i} = getUpdatedCosts(params.M{i}{params.k},model.SM);
        else
            model.KM{i} = getUpdatedCosts(params.M{i},model.SM);
        end
    end
    toc;

    %% Test the subsequence model
    model.bestThs = params.bestThs;
    model.nThreshs = params.nThreshs;
    model.scoreMeasure = params.scoreMeasure;
    model.maxWlen = params.maxWlen;
    model.k = params.k;
    display('Training the model parameters ...');
    [model,s,~] = g(model,XtrainT,Ytrain);
end