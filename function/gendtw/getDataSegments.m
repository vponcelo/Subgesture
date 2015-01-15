function [Xtrain_k,Xval_k,XLearn,I,Y,segTrain,segVal] = getDataSegments(X,Y,nSegs,k0,nmin,nmax)
% Obtain data segments
% Output:
%   Xtrain_k: training data sequences 
%   Xtrain_k: training data sequences 
%   XLearn: learning (training+validation) data sequences
%   I: Individual for the genetic algorithm
%   Y: data labels
%   segTrain: initial segmentation for training
%   segVal: initial segmentation for validation
% Input:
%   X: Data: {1} training , {2} validation , {3} test
%   Y: data labelss: {1} training , {2} validation , {3} test
%   nSegs: Number of segments to split the training sequence
%   k0: initial data clusters 
%   nmin:  minimum subsequence width
%   nmax: maximum sequence width


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
        [Xtrain,Xval] = split2GT(X,Y);
    elseif strcmp(sampling,'segments')
        [Xtrain,Xval,Y,I] = obtainSegments(X,Y,nSegs,nmin,nmax);
        I = [k0 I];
    else
        error('Specified samplign method is not valid');
    end
    Xtrain_k = [];
    Xval_k = [];
    segTrain = []; 
    segVal = [];     
end

%% Get data segments

% Coment this until %%%%%%% to sample data for each k=#gestures
if nSegs
    f = nSegs;    
else
    f = k0;
    I = [];
end
for i = 1:f
    Xtrain_k = [Xtrain_k Xtrain(Y{1}.L0 == i)];
    segTrain = [segTrain i*ones(1,length(Xtrain(Y{1}.L0 == i)))];
%     if length(Y) > 1
%         Xval_k = [Xval_k Xval(Y{2}.L0 == i)];
%         segVal = [segVal i*ones(1,length(Xval(Y{2}.L0 == i)))];
%     end
end

if ~isempty(Xval_k) && ~isempty(segVal)
    XLearn = [Xtrain_k Xval_k];    
else
    XLearn = Xtrain_k;
end