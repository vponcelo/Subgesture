function [errorX,errorTest,Ctest,accursX,accursTest] = evalCentroids(test,Z,C,labelsX,labelsTest)
% Evaluate Centroids and clusters obtained with the overall learning data

% Input:
%   test: data to test the clusters against 
%   Z: Centroid sequences
%   C: cluster indices 
%   labelsX: labels of X data
%   labelsTest: labels of test data
% Output:
%   errorX: global error over X
%   errorTest: global error over test
%   Ctest: Cluster indices over test
%   accursX: global accuracy over X
%   accursTest: global accuracy over test

cmX = confusionmat(C,labelsX);

if length(Z) > length(unique(labelsX))
    % HAY QUE PENSAR CÓMO RECALCULAR EL ERROR PARA k > #gestos
    cmX(:,length(unique(labelsX))+1:end) = [];
    error('No implementation yet for subgesturing');
elseif length(Z) < length(unique(labelsX))
    error('Number of clusters must be equal or greater than the number of classes');
end
 
%% evaluate training data
accursX = zeros(1,length(Z));
for j = 1:length(Z)
    accursX(j) = cmX(j,j)/sum(cmX(:,j));
end
meanAccur = mean(accursX);
errorX = 1-meanAccur;           

%% evaluate test data
dc = inf*ones(length(Z),length(test));
for j = 1:length(Z)
    if isempty(Z{j})
        continue;
    end
    for i = 1:length(test)
        W=dtwc(test{i},Z{j},1);
        dc(j,i) = W(end,end);
    end
end 
Ctest = zeros(1,length(test));
for i = 1:length(test)
    [~,kf] = min(dc(:,i));
    Ctest(i) = kf;
end
cmTest = confusionmat(Ctest,labelsTest);

if length(Z) > length(unique(labelsX))      % k > #gestures
    cmTest(:,length(unique(labelsTest))+1:end) = [];
end

accursTest = zeros(1,length(Z));
for j = 1:length(Z)
    accursTest(j) = cmTest(j,j)/sum(cmTest(:,j));
end
meanAccur = mean(accursTest);
errorTest = 1-meanAccur;   


