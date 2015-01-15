function showClusterResults(Ctrain,Cval,labelsT,labelsV,minErrT,minErrV,posT,posV,C2train,C2val,error2T,error2V,ZT,ZV,version)
% Draw best accuracy confusion matrix and the confusion matrix of the
% intersection of clusters giving the best accuracy
% Input:
%   Ctrain: Cluster indices for training data 
%   Cval: Cluster indices for validation data 
%   labelsT: labels for training data
%   labelsV: labels for validation data
%   minErrT: minimum error on training
%   minErrV: minimum error on validation
%   posT: execution time giving the min error on training
%   posV: execution time giving the min error on validation
%   error2T: errors of cluster intersection over training
%   error2V: errors of cluster intersection over validation
%   C2train: Intersected Cluster indices for training data
%   C2val: Intersected Cluster indices for training data
%   ZT: Chosen centroids for training data
%   ZV: Chosen centroids for validation data
%   version: string with the version of the kmeans DTW algorithm


cm = confusionmat(Ctrain,labelsT);
h=figure('visible','off');
subplot(2,2,1);imagesc(cm);
title(sprintf('Min training error: %.2f at execution %d',minErrT,posT));
subplot(2,2,2);imagesc(confusionmat(Cval,labelsV));
title(sprintf('Min validation error: %.2f at execution %d',minErrV,posV));
subplot(2,2,3);imagesc(confusionmat(C2train,labelsT));
title(sprintf('Intersection min error over training: %.2f',error2T));
subplot(2,2,4);imagesc(confusionmat(C2val,labelsV));
title(sprintf('Intersection min error over validation: %.2f',error2V));
hold on;
colormap(flipud(colormap(gray)))
hold off;
dirname = sprintf('results/chalearn2013/clustering/%s/',version);
filename = sprintf('%sk%d',dirname,length(unique(labelsT)));
saveas(h,filename,'png');
close(h)

% s=1;
% bgImage=ones(480*s,640*s)*255;
% for i = 1:length(ZV)
%     frame = round(size(ZV{i},1)/2);
%     Zdraw = [ZV{i}(frame,1:end/2); ZV{i}(frame,21:end)]'*size(bgImage,2);
%     h=showSkeleton(Zdraw*s,bgImage,1);
%     filename = sprintf('%sskels_k%d_z%d',dirname,length(ZV),i);
%     saveas(h,filename,'png');
%     close(h)
% end
