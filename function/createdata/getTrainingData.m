function [Xtrain,Ytrain] = getTrainingData(X,k,Indices,Y,l)
%% Output:
%   Xtrain: Training data for the current fold
%   Ytrain: Validation labels for the current fold (supervised)

%% Input:
%   X: Data sequences in categories
%   k: k fold
%   Indices: samples k to select
%   Y: Data Labels in categories (supervised)
%   l: gesture label

if iscell(X)
    Xtrain = X{l}(Indices{l} ~= k,:);
    if exist('Y','var')
        Ytrain = Y.Lfr(Y.Lfr==l); Ytrain = Ytrain(Indices{l} ~= k);
    else
        Ytrain = [];    % unsupervised
    end
else
    Xtrain = X(Indices ~= k,:);
    if exist('Y','var'),
        Ytrain = Y.Lfr(Indices ~= k);
    else
        Ytrain = [];    % unsupervised
    end
end