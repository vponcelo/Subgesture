function Xdiscrete = discretizeData(c,X)
% Input:
%   C: X clustered data
%   X: data to discretize

% Output
%   Xdiscrete: Discrete training data

% if ~exist('discreteData.mat','file'),
    display('Discretizing data...');
        
    X = toMat(X);
    Xdiscrete = discretizeSequence(c,X);
    
%     save('discreteData.mat','Xdiscrete');
    display('Done!');
% else
%     load('discreteData.mat');
% end
