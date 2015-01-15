function [model,s] = testGenDTWLearn(X,Y,sampling,XLearn_l,useModel,model,seqTest,GTtest)
% Validate the data sequences by means of either k-means-DTW and mean DTW models.

% output:
%   model: struct with the model parameters
%      XLearn: sub sequences generated
%      SM: Subgesture models (centroid sequences) for the k gestures
%      C: cluster indices for each sample
%      D: Dissimilarity matrix between subgesture models SM
%      M: mean sequence models
%      KM: Cost matrix of the mean sequence models in terms of SM
%      KT: Cost matrix of the test sequence in terms of SM
%      bestThs: best thresholds for the k detected gestures
%      nSegs: Number of segments to split the training sequence
%      nmin: minimum subsequence width
%      nmax: maximum sequence width
%      versions: string with the versions of the kmeans DTW algorithm
%      dist: distance metric for the k-means DTW
%      k0: initial data clusters 
%      nThreshs: Number of thresholds for testing
%   s: score (default: mean overlap for the k detected gestures)
% input:
%   X: Data: {1} training , {2} validation , {3} test
%   Y: data labelss: {1} training , {2} validation , {3} test
%   sampling: data sampling type: 'random' 'kmeansdtw'
%   useModel: flag indicating whether the validation uses the euclidean 
%       distance or the learnt model
%   seqTest: test sequence
%   GTtest: frame ground truth of the test sequence

%% Obtain data segments
[Xtrain,Xval,model.XLearn,I,Y,segTrain,segVal] = getDataSegments(X,Y,model.nSegs,model.k0,model.nmin,model.nmax);
seg0 = segTrain;

%% Temporal Clustering
model.KM = []; model.KT = [];
if useModel
    % Obtain subsets using the k-means DTW algorithm  
    kf = model.k0;                                % only test the k0 clusters
    % kf = length(unique(Y{1}.L0))-(model.k0-1);  % Force k=k0:#gestures 
    % kf = length(model.XLearn)-(model.k0-1);           % test all clusters until k=#samples

    path = 'results/chalearn2013/clustering/';

    for i = 1:length(model.versions)
        v = model.versions{i};
        if isempty(v)
            error('Version specified for the k-means DTW algorithm does not exist');
        end
        dirname = strcat(path,sampling,'/',v);
        if ~exist(dirname,'dir')
            mkdir(dirname);
        end
    %     if ~exist(sprintf('%s/kmeansdtw_%s.mat',dirname,v),'file')
            if isempty(model.XLearn)
                [CsTrain,CsVal,mErrsT,mErrsV,timeT,timeV,Z] = ...
                    runKMeansDTW(v,model.k0,model.dist,kf,Xtrain,Xval,Y,segTrain,segVal,Xtrain,Xval);
    %             save(sprintf('%s/kmeansdtw_%s.mat',dirname,v),'Xtrain','Xval','CsTrain','CsVal','segTrain','segVal','mErrsT','mErrsV','timeT','timeV','Z');
            else
                [CsTrain,~,mErrsV,~,timeV,~,Z] = runKMeansDTW(v,model.k0,model.dist,kf,[],[],Y,seg0,[],model.XLearn,[]);
    %             save(sprintf('%s/kmeansdtw_%s.mat',dirname,v),'model.XLearn','Cs','seg0','mErrsV','timeV','Z');
            end
    %     else
    %         load(sprintf('%s/kmeansdtw_%s.mat',dirname,v));
    %     end
    end

    [~,kV] = min(mErrsV);
    model.k = model.k0+kV-1;
    if isempty(model.XLearn)
        [~,kT] = min(mErrsT);
    end

    % plotClusters(dirname,model.XLearn,Xtrain,Cs,Z,Y,kf,model.k,timeV,model.k,version);

    %% Get clustered training/learning data structures     
    model.C = CsTrain{kV}{timeV(kV)};

    %% Obtain Subgesture Model for training/learning data
    if model.nSegs
        fprintf('Obtaining mean Subgesture Models from training set using k=%d ...\n',model.k);
        model.SM = Z{kV}{timeV};
        emptyCells = cellfun(@isempty,model.SM);
        model.SM(emptyCells) = [];
        display('Done!');
    end

    %% Get clustered validation data structures 
    % Only possible if clustering has been validated (so kf-model.k0+1=kV=m)
    % if CsVal{kV}{timeV(kV)}
    %     if ~model.nSegs
    %         Xval = []; labelsVal = []; 
    %         for i = 1:model.k
    %             Xval = [Xval Xval(Y{2}.L0 == i)];
    %             labelsVal = [labelsVal i*ones(1,length(Xval(Y{2}.L0 == i)))];
    %         end
    %     end
    %     if CsVal{kV}{timeV(kV)}
    %         C = CsVal{kV}{timeV(kV)};
    %     end
    %     XvalC = cell(1,model.k);
    %     for i = 1:length(XvalC)
    %         XvalC{i} = Xval(C == i);
    %     end
    %     %% Obtain Subgesture Model for validation data
    %     if model.nSegs
    %         fprintf('Computing mean Subgesture Models from validation set using k=%d ...\n',model.k);
    %         tic;            
    %         SMv = getMedianModels(XvalC,i,'direct','nogmm');
    %         toc;
    %         display('Done!');
    %     end
    % end

    %% compute Similarity matrix
    if model.nSegs
        model.D = getSimilarities(model.SM);
        if sum(~any(model.D))
            error('Some elements of the dissimilarity matrix are wrong');
        end
    end
end

%% Compute the costs of the test sequence in terms of SM 
if ~isempty(model.D)
    display('Computing the costs of the test sequence in terms of SM ...');
    model.KT = getUpdatedCosts(seqTest,model.SM);
    display('Done!');
end

%% Get median models from training/learning data
if ~any(strcmp('M',fieldnames(model)))
    if any(strcmp('k',fieldnames(model)))
        fprintf('Computing mean gesture Models from training for the m=%d distinct gestures ...\n',model.k);
        tic;
        model.M = getMedianModels(XLearn_l,model.k,'direct','nogmm');
        toc;
        display('Done!');
    else
        error('Number of clusters k has not been set');
    end
end

%% Get mean models from validation data
% Xval_l = getGroupedGestures(X,Y,2);  % it should go outside this function
% fprintf('Computing mean gesture Models from validation for the m=%d distinct gestures ...\n',model.k);
% tic;
% Mval = getMedianModels(Xval_l,model.k,'direct','nogmm');
% toc;
% display('Done!');

%% Compte costs for gesture Models M with respect to the Subgesture Models model.SM
if useModel
    if model.nSegs
        model.KM = cell(1,length(model.M));        
        for i = 1:length(model.KM)            
            model.KM{i} = getUpdatedCosts(model.M{i},model.SM);
        end
    end
end


%% Learn threshold cost parameters for each gesture
folds = 1;
% For each fold, seq (and GTtest) should be splitted in random subsets, taking
% 1 for testing
% Experiment 1: Use a learning sequence of each gesture as reference models
% (it should appear k detections with cost 0)
% Experiment 2: Use a validation sequence of each gesture as reference models 
% Experiment 3: Use Mean aligned of each gesture as reference models 
for experiment = 3:3
    fprintf('Running Experiment %d ... \n',experiment);
    tic;
    thresholds = cell(1,length(model.M));
    overlaps = cell(1,length(model.M));
    if experiment == 1        
        for k = 1:model.k
            idxsseqk = find(GTtest == k);
            i = 2; gap = false; iseqk = idxsseqk(1);
            while i < length(idxsseqk) && ~gap
                if idxsseqk(i)-idxsseqk(i-1) == 1
                      iseqk = [iseqk idxsseqk(i)];
                else
                    gap = true;
                end
                i = i + 1;
            end
            seqTest(idxsseqk(~ismember(idxsseqk,iseqk)),:) = [];
            GTtest(idxsseqk(~ismember(idxsseqk,iseqk))) = [];
        end
    end
    
    for k = 1:length(model.M)
        fprintf('Testing threshold cost parameters for gesture %d ...\n',k);
        GTtestk = GTtest == k;        
        switch experiment
            case 1        
                idxsseqk = find(GTtest == k);
                i = 2; gap = false; iseqk = idxsseqk(1);
                while i <= length(idxsseqk) && ~gap
                    if idxsseqk(i)-idxsseqk(i-1) == 1
                        iseqk = [iseqk idxsseqk(i)];
                    else
                        gap = true;
                    end
                    i = i + 1;
                end
                seqk = seqTest(iseqk,:);
                W = dtw2(seqTest,seqk,false);
            case 2
                seqVal = XLearn_l{k}{randi([1 length(XLearn_l{k})])};
                W = dtw2(seqTest,seqVal,false);
            case 3
                if useModel
                    W = dtwc(seqTest,model.M{k},false,Inf,model.D,model.KM{k},model.KT);
                    TOL_THRESH = 0.001;
                else
                    W = dtwc(seqTest,model.M{k},false);
                    TOL_THRESH = 0.01;
                end
        end
        interv = (max(W(end,2:end))-min(W(end,2:end)))/model.nThreshs;
        while interv < TOL_THRESH && model.nThreshs > 1
            model.nThreshs = round(model.nThreshs/2);
            interv = (max(W(end,2:end))-min(W(end,2:end)))/model.nThreshs;
        end
        thresholds{k} = zeros(folds,model.nThreshs);
        overlaps{k} = zeros(folds,model.nThreshs);
        for K = 1:folds        
            tMin = min(W(end,2:end));
            detSeqLog = zeros(1,length(seqTest));
            for i = 1:model.nThreshs
                thresholds{k}(K,i) = tMin + (i-1)*interv;
                idx = find(W(end,:) <= thresholds{k}(K,i));                
                for j = 1:length(idx)
                    if detSeqLog(idx(j)-1) == 0
%                         fprintf('testing with threshold %.2f\n',idx(j));
                        [in,fi,~] = aligngesture(seqTest,W(:,1:idx(j)));
                        if length(in) > 1 || length(fi) > 1
                            error('Start and end of the gesture must be scalars');
                        end
                        detSeqLog(in:fi) = 1;
                    end
                end
                overlaps{k}(K,i) = sum(GTtestk & detSeqLog)/sum(GTtestk | detSeqLog);                
            end
        end    
    end
%     display('Saving results ...');
%     filename = sprintf('%s/Exp%s.mat',outputPath,num2str(experiment));
%     save(filename,'overlaps','thresholds');
%     showOverlaps(overlaps,thresholds);    
%     display('Done!');
    toc;
end
display('Done!');

bestOverlaps = zeros(1,length(overlaps));
model.bestThs = zeros(1,length(thresholds));
for i = 1:k
    [bestOverlaps(i),pos] = max(overlaps{i});
    model.bestThs(i) = thresholds{i}(pos);
end
s = mean(bestOverlaps);
