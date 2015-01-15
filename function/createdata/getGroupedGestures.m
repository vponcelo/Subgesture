function XLearn_l = getGroupedGestures(X,Y,datapart)
%  Return a cell array with all sequences of each class in each position
% Output:
%   XLearn_l: cell with grouped sequences
% Input:
%   X: data sequences
%   Y: data labels
%   datapart: partition of the data to consider: 0: all , 1: only training,
%       2: only validation

if length(Y) > 1
    switch datapart
        case 0
            labelsLearn = [Y{1}.L Y{2}.L];
        case 1
            labelsLearn = [Y{1}.L];
        case 2
            labelsLearn = [Y{2}.L];
    end    
elseif length(Y) == 1
    display('Only training data partition can be considered');
    labelsLearn = Y{1}.L;
else
    error('There are no data labels');
end

XLearn_l = cell(1,length(unique(labelsLearn)));
for i = 1:length(XLearn_l)
    XLearn_l{i} = cell(1,sum(labelsLearn==i)); 
    idxLval = find(labelsLearn == i);
    k = datapart;
    if datapart == 0, k = 1; end
    for j = 1:length(XLearn_l{i})
        if idxLval(j) > length(Y{k}.L)
            idxLval(j:end) = idxLval(j:end) - length(Y{k}.L);
            k = k + 1;
        end
        XLearn_l{i}{j} = X{k}(Y{k}.seg(idxLval(j)):Y{k}.seg(idxLval(j)+1)-1,:);        
    end
end