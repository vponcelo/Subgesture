function [Ctrain,Cval,errT,errV] = intersClust(CsTrain,CsVal,accursT,accursV,labelsT,labelsV)
% Compute the intersection of clusters giving the best accuracy

% Input:
%   CsTrain: Cluster indices for training data
%   CsVal: Cluster indices for validation data
%   accursT: all accuracies over training
%   accursV: all accuracies over validation
%   labelsT: labels of training data
%   labelsV: labels of validation data
% Output:
%   Ctrain: Intersected cluster indices for training
%   Cval: Intersected cluster indices for validation
%   errT: Intersection min error for training
%   errV: Intersection min error for validation

%% Compute cluster intersection over training 
bestCs = zeros(size(accursT,2),2);
Ctrain = [];
for j = 1:size(bestCs,1)
    [bestCs(j,1),bestCs(j,2)] = max(accursT(:,j));    
    idx = CsTrain{bestCs(j,2)} == j;
    if isempty(Ctrain)
        Ctrain = j*idx;
    else
        for i = 1:length(Ctrain)
            if Ctrain(i) > 0 && j*idx(i) == j
                Ctrain(i) = j*idx(i);
            elseif Ctrain(i) == 0 
                Ctrain(i) = j;
            end
        end
    end
end
cm = confusionmat(Ctrain,labelsT);
bestAccur = sum(diag(cm))/sum(sum(cm));
errT = 1-bestAccur;

%% Compute cluster intersection over validation
bestCs = zeros(size(accursV,2),2);
Cval = [];
for j = 1:size(bestCs,1)            
    [bestCs(j,1),bestCs(j,2)] = max(accursV(:,j));    
    idx = CsVal{bestCs(j,2)} == j;
    if isempty(Cval)
        Cval = j*idx;
    else
        for i = 1:length(Cval)
            if Cval(i) > 0 && j*idx(i) == j
                Cval(i) = j*idx(i);
            elseif Cval(i) == 0 
                Cval(i) = j;
            end
        end
    end
end
cm = confusionmat(Cval,labelsV);
bestAccur = sum(diag(cm))/sum(sum(cm));
errV = 1-bestAccur;


