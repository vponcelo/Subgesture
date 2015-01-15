function X = norm2neck(X)

global COORDS;
global JOINTS;

if strcmp(COORDS,'world')
    oc = 3;
elseif strcmp(COORDS,'pixel')
    oc = 2;
end

[head,pos] = ismember(4,JOINTS);
if ~head
    error('norm2neck:normError','normalization is not possible because neck joint is not contained in the data');
end
if iscell(X)
    for i = 1:length(X)
        X{i}(:,1:length(JOINTS)) = X{i}(:,1:length(JOINTS)) - repmat(X{i}(:,pos),1,length(JOINTS));
        X{i}(:,length(JOINTS)+1:length(JOINTS)*2) = X{i}(:,length(JOINTS)+1:length(JOINTS)*2) - repmat(X{i}(:,length(JOINTS)+pos),1,length(JOINTS));
        if strcmp(COORDS,'world')
            X{i}(:,length(JOINTS)*2+1:length(JOINTS)*3) = X{i}(:,length(JOINTS)*2+1:length(JOINTS)*3) - repmat(X{i}(:,length(JOINTS*2)+pos),1,length(JOINTS));
        end
%         for j = 1:size(X{i},1)
%             X{i}(:,1:length(JOINTS)) = X{i}(:,1:length(JOINTS)) - X{i}(:,pos);
%             X{i}(:,length(JOINTS)+1:length(JOINTS)*2) = X{i}(:,length(JOINTS)+1:length(JOINTS)*2) - X{i}(:,length(JOINTS)+pos);
%             if strcmp(COORDS,'world')
%                 X{i}(:,length(JOINTS)*2+1:length(JOINTS)*3) = X{i}(:,length(JOINTS)*2+1:length(JOINTS)*3) - X{i}(:,length(JOINTS)*2+pos);
%             end            
%         end
        %% Data normalization by each position feature
        minValues = repmat(min(X{i}),size(X{i},1),1);
        maxValues = repmat(max(X{i}),size(X{i},1),1);
        X{i} = (X{i}-minValues) ./ (maxValues-minValues); 
        X{i}(:,pos:length(JOINTS):length(JOINTS)*oc) = 0;
        X{i}(:,any(isnan(X{i})))=0;
    end
end