function XLearn_l = getGroupedGestures(X,Y,datapart,Xtest,Ytest)
%  Return a cell array with all sequences of each class in each position
% Output:
%   XLearn_l: cell with grouped sequences
% Input:
%   X: data sequences
%   Y: data labels
%   datapart: partition of the data to consider: 0: all , 1: only training,
%       2: only validation

display('Obtaining gestures grouped by gesture class...');

if length(Y) > 1 || ~isempty(Ytest)
    switch datapart
        case 0
            labelsLearn = [Y{1}.L Y{2}.L];
        case 1
            labelsLearn = [Y{1}.L];
        case 2
            labelsLearn = [Y{2}.L];
        case 3
            labelsLearn = [Ytest.L];
    end    
elseif length(Y) == 1
    display('Only training data partition can be considered');
    labelsLearn = Y{1}.L;
else
    error('There are no data labels');
end

XLearn_l = cell(1,length(unique(labelsLearn)));
for i = 1:length(XLearn_l)
    k = datapart;
    if datapart == 0, k = 1; end
    if k == 3, Yc = Ytest; Xc = Xtest; elseif k < 3, Yc = Y{k}; Xc = X{k}; end
    if length(fieldnames(Yc)) == 3
        XLearn_l{i} = cell(1,sum(labelsLearn==i)); 
    end
    idxLval = find(labelsLearn == i);
    if iscell(XLearn_l{i})
        for j = 1:length(idxLval)
            if idxLval(j) > length(Yc.L)
                idxLval(j:end) = idxLval(j:end) - length(Yc.L);
                k = k + 1;
            end
            XLearn_l{i}{j} = Xc(Yc.seg(idxLval(j)):Yc.seg(idxLval(j)+1)-1,:);                
        end
    else
        XLearn_l{i} = Xc(Yc.Lfr==i,:);
    end
end
display('Done!');