function [cluster,dist]=getCluster(centers,vector)
     distMat=pdist2(centers,vector);
     [dist,cluster]=min(distMat);