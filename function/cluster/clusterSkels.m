function [Idx,C] = clusterSkels( features, k, seq, skels, nframesSeg )
% Cluster the skeletons using kmeans
% 
% Input -
%   features: cell array with the skeleton sequences
%   k: number of clusters to group and visualize skeletons
%   seq: number of the sample sequence (batch number)    
%   skels: skeleton sequences 
%   nframesSeg: fixed number of frames to segment the sequences

% Output -
%   Idx: Clusters of each sample
%   C: Centroids of each k-mean cluster

    if k
        % Apply K-Means
        if iscell(features),
            features = toMat(features);
        end
        
        [Idx,C] = kmeans(features,k);        
        
        if exist('skels','var') && exist('seq','var'),
            path = 'results/chalearn2012/visualClusters/kmeans';
            % Store the skeletons
            s=0.5;
            bgImage=ones(480*s,640*s)*255;
            for c=1:k,
                cPath=sprintf('%s/%03d',path,c);
                if ~exist(cPath,'dir'),
                    mkdir(cPath);
                end
                % show skeletons falling in the cluster c
                Indices = find(Idx == c);
                for j=1:length(Indices)
                    iSeq = 1;
                    iSkel = Indices(j);
                    if nframesSeg > 1,
                        iSkel = iSkel*nframesSeg - round(nframesSeg/2);
                    end
                    while iSkel > length(skels{iSeq})
                        iSkel = iSkel - length(skels{iSeq});
                        iSeq = iSeq + 1;                        
                    end                             
                    showSkeleton(skels{iSeq}(iSkel).PixelPosition*s,bgImage,1);
                    filepath = sprintf('%s/%05d_%04d.jpg',cPath,iSeq-1+seq,iSkel);
                    if ~exist(filepath,'file'),
                        img=getframe;
                        imwrite(img.cdata,filepath); 
                        close all;                             
                    end                    
                end
            end
        end
    end
end

