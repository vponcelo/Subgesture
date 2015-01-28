function [Xtrain,Ytrain] = getTrainingData(X,k,Indices,Y)
%% Output:
%   Xtrain: Training data for the current fold
%   Ytrain: Validation labels for the current fold (supervised)

%% Input:
%   X: Data sequences in categories
%   Y: Data Labels in categories (supervised)

Xtrain = X(Indices ~= k,:);
if exist('Y','var'),
    Ytrain = Y.Lfr(Indices ~= k);
else
    Ytrain = [];
end