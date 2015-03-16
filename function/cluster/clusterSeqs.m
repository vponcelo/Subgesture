function [Z,CsX,CsVal,timeT,timeV,mErrT,mErrV] = clusterSeqs(X,Xval,k,dist,labelsT,labelsV,version,resize)
% Obtain k sequence clusters from X

% output:
%   CsX: cluster indices for all nTimes executions for training data
%   CsVal: cluster indices for all nTimes executions for validation data
%   mErrsT: mean error of training for all nTimes executions
%   mErrsV: mean error of validation for all nTimes executions
%   timeT: execution time giving the min error on training
%   timeV: execution time giving the min error on validation
%   Z: Chosen centroids for all nTimes executions
% input:
%   X: input sequence data (training)
%   Xval: input sequence data (validation)
%   k: number of clusters
%   dist: distance metric for the k-means DTW
%   labelsT: labels of training data
%   labelsV: labels of validation data
%   version: string with the version of the kmeans DTW algorithm
%   resize: Use resizing instead of mean DTW alignment
   
if strcmp(version(4),'0')
    nTimes = 1;
else
    nTimes = 1;
end

errorsT = zeros(1,nTimes);
errorsV = zeros(1,nTimes);
CsX = cell(1,nTimes);
CsVal = cell(1,nTimes);
accursT = zeros(nTimes,k);
accursV = zeros(nTimes,k);
Z = cell(1,nTimes);
% matlabpool open;
% parfor when nTimes > 1
for i = 1:nTimes
    CsX{i} = zeros(1,length(X));
    CsVal{i} = zeros(1,length(Xval));
%     fprintf(sprintf('Execution %d of the k-means DTW algorithm ...\n',i));
    [Z{i},CsX{i}] = kmeansDTW(X,k,version,dist,labelsT,resize);
    if ~isempty(Xval)
        if k == length(labelsT.L)
            [errorsT(i),errorsV(i),CsVal{i},accursT(i,:),accursV(i,:)] = ...
                evalCentroids(Xval,Z{i},CsX{i},labelsT,labelsV);
        end
    end
end
% matlabpool close;

if ~isempty(Xval)
    if k == length(labelsT.L) 
        [minErrT,timeT] = min(errorsT);

        [C2train,C2val,error2T,error2V] = intersClust(CsX,CsVal,accursT,accursV,labelsT,labelsV);
        fprintf(sprintf('Data clusters successfully obtained.\n'));

        Ctrain = CsX{timeT};
        ZT = Z{timeT};
        [minErrV,timeV] = min(errorsV);
        Cval = CsVal{timeV};
        ZV = Z{timeV};

        showClusterResults(Ctrain,Cval,labelsT,labelsV,minErrT,minErrV,timeT,...
            timeV,C2train,C2val,error2T,error2V,ZT,ZV,version);
        mErrT = mean(errorsT);
        mErrV = mean(errorsV);
    end
else
    timeT = 1; timeV = 1; mErrT = Inf; mErrV = Inf;
end
 

