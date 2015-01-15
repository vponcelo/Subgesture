function [idxRowsDel,X] = cleanData(X)

idxDel = cell(1,2);
idxRowsDel = cell(1,2);
for i = 1:length(X)
    idxDel{i} = X{i}(:,1:20) == 0 | X{i}(:,21:40) == 0;
%     sum(sum(idxDel{i}'))
    idxRowsDel{i} = all(idxDel{i}');    
end
if sum(cell2mat(idxRowsDel)) == 0
    return;
end



%% Interpolation (not fully implemented)
% for i = 1:length(X)
% %     idx = find(idxDel{i} == 1);
%     iniInterp = idx(1);
%     fiInterp = -1;
%     for j = 2:size(idxDel{i},1)
%         if idx(j)-idx(j-1) > 1
%             fiInterp = idx(j-1);
%             fiX = X{i}(fiInterp+1,:);
%             inX = X{i}(iniInterp,:);
%             yi = interp1q(X{i}(iniInterp:fiInterp,:),fiX,fiInterp-iniInterp+1)            
%         else
%             if fiInterp > 0
%                 iniInterp = idx(j);
%                 fiInterp = -1;
%             end            
%         end
%     end
% end
