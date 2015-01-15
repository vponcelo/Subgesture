function [m_overlap,m_ths] = testGenDTW(X,Y,sampling,nSegs,nmin,nmax,versions,dist,k0,valEuDist)
% Validate the data sequences by means of either k-means-DTW and mean DTW models.

% output:
%   m_overlap: mean overlap for the k detected gestures
%   m_ths: mean threshold for the k detected gestures
% input:
%   X: Data: {1} training , {2} validation , {3} test
%   Y: data labels: {1} training , {2} validation , {3} test
%   sampling: data sampling type: 'random' 'kmeansdtw'
%   nSegs: Number of segments to split the training sequence
%   nmin:  minimum subsequence width
%   nmax: maximum sequence width
%   versions: string with the versions of the kmeans DTW algorithm
%   dist: distance metric for the k-means DTW
%   k0: initial data clusters 
%   valEuDist: flag indicating whether the validation uses the euclidean
%       distance or not


%% Initial data splitting
if nSegs
    sampling = 'segments';
end
if strcmp(sampling,'random')
    % Obtain k0 random subsets
    Xtrain = getRandomSubsets(X{1},k0);
    Xval = getRandomSubsets(X{2},k0);     
else
    % Obtain gesture subsets
    if strcmp(sampling,'labels')
        Y{1}.L0 = Y{1}.L;
        Y{2}.L0 = Y{2}.L;
        Y{1}.seg0 = Y{1}.seg;
        Y{2}.seg0 = Y{2}.seg;
        [Xtrain,Xval] = splitIntoSegs(X,Y);        
    elseif strcmp(sampling,'segments')
        [Xtrain,Xval,Y] = obtainSegments(X,Y,nSegs,nmin,nmax);               
    else
        error('Specified samplign method is not valid');
    end
    Xtrain_k = [];
    Xval_k = [];
    labelsTrain = []; 
    labelsVal = [];     
end

%% Get data grouped by gestures
% This should be the normal way (outside the loop), where data should not change.

% Coment this until %%%%%%% to sample data for each k=#gestures
if nSegs
    f = nSegs;
else
    f = k0;
end
for i = 1:f
    Xtrain_k = [Xtrain_k Xtrain(Y{1}.L0 == i)];
    labelsTrain = [labelsTrain i*ones(1,length(Xtrain(Y{1}.L0 == i)))];
end
kf = length(Xtrain_k)-(k0-1);     % test all clusters until k=#samples
for i = 1:f
    Xval_k = [Xval_k Xval(Y{2}.L0 == i)];
    labelsVal = [labelsVal i*ones(1,length(Xval(Y{2}.L0 == i)))];
end 
%%%%%%%

kV = 1;

%% Temporal Clustering
% Obtain subsets using the k-means DTW algorithm  
%   Comment lines 47-54 to obtain data subsets at each k until 
%   k=#gestures and uncomment line 65 or 66 as desired
kf = k0;                               % only test 1 cluster k0 
% kf = length(unique(Y{1}.L0))-(k0-1);    % Force k=k0:#gestures 
    
path = 'results/chalearn2013/clustering/';
    
for i = 1:length(versions)
    v = versions{i};
    if isempty(v)
        error('Version specified for the k-means DTW algorithm does not exist');
    end
    dirname = strcat(path,sampling,'/',v);
    if ~exist(dirname,'dir')
        mkdir(dirname);
    end
%     if ~exist(sprintf('%s/kmeansdtw_%s.mat',dirname,v),'file')
        [CsTrain,CsVal,mErrsT,mErrsV,timeT,timeV,Z] = ...
            runKMeansDTW(v,k0,dist,kf,Xtrain,Xval,Y,labelsTrain,labelsVal,Xtrain_k,Xval_k);
%         save(sprintf('%s/kmeansdtw_%s.mat',dirname,v),'Xtrain_k','Xval_k','CsTrain','CsVal','labelsTrain','labelsVal','mErrsT','mErrsV','timeT','timeV','Z');
%     else
%         load(sprintf('%s/kmeansdtw_%s.mat',dirname,v));
%     end
end

[~,kT] = min(mErrsT);
[~,kV] = min(mErrsV);

% for i = kf-k0+1
%     createSeqVideo(Z{i}{timeV(i)},dirname);
% end
% for i = 1:kf-k0+1
%     outputPath = sprintf('%s/C%d',dirname,i);
%     if ~exist(outputPath,'dir');
%         mkdir(outputPath);
%     end
%     createSeqVideo(Xtrain(CsTrain{kV}{timeV(kV)}==i),outputPath);        
% end
% h=figure('visible','off');
% subplot(2,1,1);plot(k0:length(unique(Y{1}.L0)),mErrsT);
% title(sprintf('Version %s: Mean clustering errors over training data',versions));
% subplot(2,1,2);plot(k0:length(unique(Y{1}.L0)),mErrsV);
% title(sprintf('Version %s: Mean clustering errors over validation data',versions));
% hold on
% filename = sprintf('%s/clusterErrors.png',dirname);
% if ~exist(filename,'file')
%     saveas(h,filename,'png');
% end
% hold off;
% close(h);

%% Get clustered training data structures     
if ~nSegs
    Xtrain_k = []; labelsTrain = [];
    for i = 1:k0+kV-1
        Xtrain_k = [Xtrain_k Xtrain(Y{1}.L0 == i)];
        labelsTrain = [labelsTrain i*ones(1,length(Xtrain(Y{1}.L0 == i)))];
    end
end
C = CsTrain{kV}{timeV(kV)};

%% Obtain Subgesture Model for training data
if nSegs
    fprintf('Obtaining mean Subgesture Models from training set using k=%d ...\n',k0+kV-1);
    SMt = Z{kV}{timeV};
    display('Done!');
end

%% Get clustered validation data structures 
% Only possible if clustering has been validated (so kf-k0+1=kV=m)
% if CsVal{kV}{timeV(kV)}
%     if ~nSegs
%         Xval_k = []; labelsVal = []; 
%         for i = 1:k0+kV-1
%             Xval_k = [Xval_k Xval(Y{2}.L0 == i)];
%             labelsVal = [labelsVal i*ones(1,length(Xval(Y{2}.L0 == i)))];
%         end
%     end
%     if CsVal{kV}{timeV(kV)}
%         C = CsVal{kV}{timeV(kV)};
%     end
%     XvalC = cell(1,k0+kV-1);
%     for i = 1:length(XvalC)
%         XvalC{i} = Xval_k(C == i);
%     end
%     %% Obtain Subgesture Model for validation data
%     if nSegs
%         fprintf('Computing mean Subgesture Models from validation set using k=%d ...\n',k0+kV-1);
%         tic;            
%         SMv = getMedianModels(XvalC,i,'direct','nogmm');
%         toc;
%         display('Done!');
%     end
% end

%% compute Similarity matrix
if nSegs
    if valEuDist
        D = getSimilarities(SMt);
    else
        D = [];
    end
end

%% Get median models from training data
% Xtrain_l = cell(1,length(unique(Y{1}.L)));
% for i = 1:length(unique(Y{1}.L))
%     Xtrain_l{i} = cell(1,sum(Y{1}.L == i)); 
%     idxLval = find(Y{1}.L == i);    
%     for j = 1:length(Xtrain_l{i})
%         Xtrain_l{i}{j} = X{1}(Y{1}.seg(idxLval(j)):Y{1}.seg(idxLval(j)+1)-1,:);
%     end
% end
% fprintf('Computing mean gesture Models from training for the m=%d distinct gestures ...\n',length(unique(Y{1}.L)));
% tic;
% Mtrain = getMedianModels(Xtrain_l,i,'direct','nogmm');
% toc;
% display('Done!');

%% Get mean models from validation data
Xval_l = cell(1,length(unique(Y{2}.L)));
for i = 1:length(unique(Y{2}.L))
    Xval_l{i} = cell(1,sum(Y{2}.L == i)); 
    idxLval = find(Y{2}.L == i);    
    for j = 1:length(Xval_l{i})
        Xval_l{i}{j} = X{2}(Y{2}.seg(idxLval(j)):Y{2}.seg(idxLval(j)+1)-1,:);
    end
end
fprintf('Computing mean gesture Models from validation for the m=%d distinct gestures ...\n',k0-kV+1);
tic;
Mval = getMedianModels(Xval_l,k0-kV+1,'direct','nogmm');
toc;
display('Done!');

%% Compte costs for gesture Models M with respect to the Subgesture Models SM
if nSegs
    KM = cell(1,length(Mval));
    for i = 1:length(KM)
        KM{i} = getUpdatedCosts(Mval{i},SMt);
    end
end

%% Validation 
% Get a sequence to test against
fps = 20;
secs = 60;
l = [24 78 150]; % 78 more samples for each gesture when k=3
GTseq = cell(1,length(l));
seq = cell(1,length(l));
for i = 1:length(l)
    j = l(i);
    GTseq{i} = zeros(1,fps*secs);
    seq{i} = zeros(fps*secs,size(X{1},2));
    % Obtain continuous sequence of time fps*secs (aprox)
    while Y{1}.seg(j) < Y{1}.seg(l(i))+fps*secs-1
        j = j + 1;
        if j > l(i)+1
            startSeq = endSeq + 1;
        else
            startSeq = 1;
        end    
        endSeq = startSeq + (Y{1}.seg(j)-1-Y{1}.seg(j-1));
        seq{i}(startSeq:endSeq,:) = X{1}(Y{1}.seg(j-1):Y{1}.seg(j)-1,:);
        GTseq{i}(startSeq:endSeq) = Y{1}.L(j-1);    
    end
end
GTseq = cell2mat(GTseq); seq = toMat(seq);
idxDel = GTseq >= 13 | GTseq == 0;
GTseq(idxDel) = []; seq(idxDel,:) = [];
% outputPath = 'results/chalearn2013/validation';
% createSeqVideo(seq,outputPath);      % create video sequence

%% Compute costs for test sequence with respect to the Subgesture Models SM
KT = getUpdatedCosts(seq,SMt);

%% Learn threshold cost parameters for each gesture
folds = 1;
% For each fold, seq (and GTseq) should be splitted in random subsets, taking
% 1 for testing
nThresholds = 100;
% Experiment 1: Use a training sequence of each gesture as reference models
% (it should appear k detections with cost 0)
% Experiment 2: Use a validation sequence of each gesture as reference models 
% Experiment 3: Use Mean aligned of each gesture as reference models 
for experiment = 3:3
    fprintf('Running Experiment %d ... \n',experiment);
    tic;
    thresholds = cell(1,k0);
    overlaps = cell(1,k0);
    if experiment == 1        
        for k = 1:k0+kV-1
            idxsseqk = find(GTseq == k);
            i = 2; gap = false; iseqk = idxsseqk(1);
            while i < length(idxsseqk) && ~gap
                if idxsseqk(i)-idxsseqk(i-1) == 1
                      iseqk = [iseqk idxsseqk(i)];
                else
                    gap = true;
                end
                i = i + 1;
            end
            seq(idxsseqk(~ismember(idxsseqk,iseqk)),:) = [];
            GTseq(idxsseqk(~ismember(idxsseqk,iseqk))) = [];
        end
    end
    
    for k = 1:k0+kV-1
        thresholds{k} = zeros(folds,nThresholds);
        overlaps{k} = zeros(folds,nThresholds);
        fprintf('Learning threshold cost parameters for gesture %d ...\n',k);
        GTseqk = GTseq == k;        
        switch experiment
            case 1        
                idxsseqk = find(GTseq == k);
                i = 2; gap = false; iseqk = idxsseqk(1);
                while i <= length(idxsseqk) && ~gap
                    if idxsseqk(i)-idxsseqk(i-1) == 1
                        iseqk = [iseqk idxsseqk(i)];
                    else
                        gap = true;
                    end
                    i = i + 1;
                end
                seqTrain = seq(iseqk,:);
                W = dtw2(seq,seqTrain,false);
            case 2
                seqVal = Xval_l{k}{randi([1 length(Xval_l{k})])};
                W = dtw2(seq,seqVal,false);
            case 3
                W = dtw2(seq,Mval{k},false,Inf,D,KM{k},KT);
        end
        interv = (max(W(end,2:end))-min(W(end,2:end)))/nThresholds;
        for K = 1:folds        
            fCost = 0;
            tMin = min(W(end,2:end));
            detSeqLog = zeros(1,length(seq));
            for i = 1:nThresholds
                thresholds{k}(K,i) = tMin + (i-1)*interv;
                idx = find(W(end,:) <= thresholds{k}(K,i));
                idx = idx(idx > fCost);
                for j = 1:length(idx)
%                     fprintf('testing with threshold %.2f\n',idx(j));
                    % idx(1): primer indice de la matriz W donde hay coste < c
                    [in,fi,~] = aligngesture(seq,W(:,1:idx(j)));
                    if length(in) > 1 || length(fi) > 1
                        error('Start and end of the gesture must be scalars');
                    end
                    detSeqLog(in:fi) = 1;              
                end
                overlaps{k}(K,i) = sum(GTseqk & detSeqLog)/sum(GTseqk | detSeqLog);
                if ~isempty(idx)
                    fCost = idx(end);
                end                
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

bestOverlaps = zeros(1,k);
bestThs = zeros(1,k);
for i = 1:k
    [bestOverlaps(i),pos] = max(overlaps{i});
    bestThs(i) = thresholds{i}(pos);
end
m_overlap = mean(bestOverlaps);
m_ths = mean(bestThs);