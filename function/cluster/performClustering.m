function Ctrain = performClustering(X,Y,type,numClusters,numIterations)
% Input:
%   X: Training data
%   Y: Training data labels (empty list when unsupervised)
%   type: clustering method
%   numClusters: number of centroids to clusterize the data
%   numIterations: number of iterations for the clustering process

global DATATYPE;

% Output:
%   Ctrain: training data clusters

if ~exist(strcat('results/',DATATYPE,'/clustering/clusters.mat'),'file'),
    display('Clustering data...');

    if numClusters == 0,
        len = 0;
        for i = 1:length(X),
            len = len + length(X{i});
        end
        numClusters = round(len/100);
    end
    
    %% Clustering methods for SubGesture categorization
    switch type
        case 'kmlsample'           
            %% KMLSample clustering
            display('Method: kmlsample (k-means based)');
            if ~exist('numIterations','var'),
                error('Max number of iterations is needed for the KML clustering method');
            end
            Ctrain = getDataClusters(numClusters,numIterations,X);
        case 'kmeans'
            %% k-means clustering
            display('Method: Mathworks k-means');
            [~,Ctrain] = clusterSkels( X, numClusters );     
        case 'haca'
            %% Hierarchical Aligned Clustering Analysis
            doHaca(X,Y);
    end
    
    save(strcat('results/',DATATYPE,'/clustering/clusters.mat'),'Ctrain');

    display('Done!');
else
    load(strcat('results/',DATATYPE,'/clustering/clusters.mat'));
end