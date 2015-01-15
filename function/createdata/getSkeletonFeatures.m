function [features,skeletons] = getSkeletonFeatures( skeletons )
    % Input -
    %  allSkeletons: cell array with the skeleton sequences
    
    % Output -
    %   features: features of concatenated skeletons 
    
    global COORDS;
    global JOINTS;
    
    features = cell(1,length(skeletons));
    for k = 1:length(skeletons)
        features{k} = cell(1,length(skeletons{k}));    
        for i = 1:length(features{k})
            nskels = length(skeletons{k}{i});
            if nskels == 1
                nskels = length(skeletons{k}{i}.skeleton);
            end
            if strcmp(COORDS,'world')
                features{k}{i} = zeros(nskels,length(JOINTS)*3);
            elseif strcmp(COORDS,'pixel')
                features{k}{i} = zeros(nskels,length(JOINTS)*2);
            end
            for j=1:nskels,
                if length(skeletons{k}{i}) == 1
                    if strcmp(COORDS,'world')
%                         features{k}{i}(j,:)=reshape(skeletons{k}{i}.skeleton(j).WorldPosition,1,size(features{k}{i},2));
                        features{k}{i}(j,:)=reshape(skeletons{k}{i}.skeleton(j).WorldPosition(JOINTS,:),1,size(features{k}{i},2));
                        features2{k}{i}(j,:)=skeletons{k}{i}.angles(j).BufH;
                    elseif strcmp(COORDS,'pixel')
%                         features{k}{i}(j,:)=reshape(skeletons{k}{i}.skeleton(j).PixelPosition,1,size(features{k}{i},2));
                        features{k}{i}(j,:)=reshape(skeletons{k}{i}.skeleton(j).PixelPosition(JOINTS,:),1,size(features{k}{i},2));
                        features2{k}{i}(j,:)=skeletons{k}{i}.angles(j).BufH;
                    end                    
                else
                    if strcmp(COORDS,'world')
%                         features{k}{i}(j,:)=reshape(skeletons{k}{i}(j).WorldPosition,1,size(features{k}{i},2));
                        features{k}{i}(j,:)=reshape(skeletons{k}{i}(j).WorldPosition(JOINTS,:),1,size(features{k}{i},2));
                        features2{k}{i}(j,:)=skeletons{k}{i}.angles(j).BufH;
                    elseif strcmp(COORDS,'pixel')
%                         features{k}{i}(j,:)=reshape(skeletons{k}{i}(j).PixelPosition,1,size(features{k}{i},2));
                        features{k}{i}(j,:)=reshape(skeletons{k}{i}(j).PixelPosition(JOINTS,:),1,size(features{k}{i},2));
                        features2{k}{i}(j,:)=skeletons{k}{i}.angles(j).BufH;
                    end
                end
            end
            %features{k}{i} = [features{k}{i} skeletons{k}{i}.skeleton(j).WorldRotation];
%             Indices = features{k}{i} == 0;    
%             features{k}{i}(Indices(:,1)==1,:) = [];
%             if length(skeletons{k}{i}) == 1
%                 skeletons{k}{i}.skeleton(Indices(:,1)==1) = [];
%             else
%                 skeletons{k}{i}(Indices(:,1)==1) = [];
%             end
        end
%         features{k} = toMat(features{k});
        features{k} = [toMat(features{k}),toMat(features2{k})];
    end
end
