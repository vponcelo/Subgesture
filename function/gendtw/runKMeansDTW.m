function [CsTrain,CsVal,mErrsT,mErrsV,timeT,timeV,Z] = ...
    runKMeansDTW(params,Xtrain,Xval,Y,labelsT,labelsV,Xtrain_k,Xval_k)
% Prepare the data and parameters to run the specified version of the  
% k-means-DTW algorithm
% Input:
%   params: set of parameters
%   Xtrain: Training data 
%   Xval: Validation data
%   Y: Learning data labels structure
%   labelsT: labels of training data
%   labelsV: labels of validation data
%   Xtrain_k: Training data to cluster in k splits (it is created if does not exist)
%   Xtrain_k: Validation data to cluster in k splits (it is created if does not exist)

%% Split and visualize learning data by means of k-means dtw clustering  
if isempty(Xtrain_k) && isempty(labelsT) && isempty(Xval_k) && isempty(labelsV)
    sampleData = true;
else
    sampleData = false;
end

k0=params.k0;kf=k0;v=params.version;r=params.resize;d=params.dist;
CsTrain = cell(1,kf-k0+1);
CsVal = cell(1,kf-k0+1);
mErrsT = zeros(1,kf-k0+1);
mErrsV = zeros(1,kf-k0+1);
timeT = zeros(1,kf-k0+1);
timeV = zeros(1,kf-k0+1);
Z = cell(1,kf-k0+1);

for k = k0:1:kf
    if sampleData
        Xtrain_k = []; Xval_k = []; labelsT  = []; labelsV = [];
        for i = 1:k
            Xtrain_k = [Xtrain_k Xtrain(Y{1}.L == i)];
            labelsT = [labelsT i*ones(1,length(Xtrain(Y{1}.L == i)))];        
            Xval_k = [Xval_k Xval(Y{2}.L == i)];
            labelsV = [labelsV i*ones(1,length(Xval(Y{2}.L == i)))];
        end
    end

	display(sprintf(sprintf('Obtaining data clusters using k=%d ...',k)));
%     tic;
    if ~isempty(Xval_k)
        [Z{k-k0+1},CsTrain{k-k0+1},CsVal{k-k0+1},timeT(k-k0+1),...
            timeV(k-k0+1),mErrsT(k-k0+1),mErrsV(k-k0+1)] = ...
            clusterSeqs(Xtrain_k,Xval_k,k,d,labelsT,labelsV,v,r);
    else
        [Z{k-k0+1},CsTrain{k-k0+1},~,timeT(k-k0+1),...
            ~,mErrsT(k-k0+1),~] = ...
            clusterSeqs(Xtrain_k,[],k,d,labelsT,[],v,r);
    end
%     toc;    
end

