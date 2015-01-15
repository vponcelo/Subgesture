function [X_train,X_val,Y_train,Y_val] = getLearningData(X,k,Indices,Y)
%% Output:
%   X_train: Training data for the current fold
%   X_train: Validation data for the current fold
%   X_val: Validation data for the current fold
%   Y_val: Validation labels for the current fold (supervised)

%% Input:
%   X: Data sequences in categories
%   Y: Data Labels in categories (supervised)

X_val = X(Indices == k);
X_train = X(Indices ~= k);
Y_train = [];
Y_val = [];
if exist('Y','var'),
    Y_val = Y(Indices == k);
    Y_train = Y(Indices ~= k);
else
    Y_val = [];
    Y_train = [];
end