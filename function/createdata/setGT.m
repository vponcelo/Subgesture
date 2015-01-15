function Y = setGT(X,Y,nSegs,nframesSeg)
% Add / Create the label vector to Y with the Ground Truth

% Output -
%   Y: Label cell containing the label vector for each sample
% Input - 
%   X: Input data
%   Y: Data segments and labels
%   nSegs: number of segments
%   nframesSeg: fixed number of frames to segment the sequences '0' means no fixed segmentation ; '1' means frame segmentation

if ~isempty(Y)
    cnames = Y.cnames;
    Y.L = zeros(1,length(cnames));
    i = 1;
    usedLabels = cnames;
    while ~isempty(usedLabels) 
        labels = strcmp(usedLabels{1}, cnames);
        if strcmp(usedLabels{1}, 'vattene');
            Y.L(labels) = 1;       
        elseif strcmp(usedLabels{1}, 'vieniqui');
            Y.L(labels) = 2;       
        elseif strcmp(usedLabels{1}, 'perfetto');
            Y.L(labels) = 3;       
        elseif strcmp(usedLabels{1}, 'furbo');
            Y.L(labels) = 4;       
        elseif strcmp(usedLabels{1}, 'cheduepalle');
            Y.L(labels) = 5;       
        elseif strcmp(usedLabels{1}, 'chevuoi');
            Y.L(labels) = 6;       
        elseif strcmp(usedLabels{1}, 'daccordo');
            Y.L(labels) = 7;       
        elseif strcmp(usedLabels{1}, 'seipazzo');
            Y.L(labels) = 8;       
        elseif strcmp(usedLabels{1}, 'combinato');
            Y.L(labels) = 9;       
        elseif strcmp(usedLabels{1}, 'freganiente');
            Y.L(labels) = 10;       
        elseif strcmp(usedLabels{1}, 'ok');
            Y.L(labels) = 11;       
        elseif strcmp(usedLabels{1}, 'cosatifarei');
            Y.L(labels) = 12;       
        elseif strcmp(usedLabels{1}, 'basta');
            Y.L(labels) = 13;       
        elseif strcmp(usedLabels{1}, 'prendere');
            Y.L(labels) = 14;       
        elseif strcmp(usedLabels{1}, 'noncenepiu');
            Y.L(labels) = 15;       
        elseif strcmp(usedLabels{1}, 'fame');
            Y.L(labels) = 16;       
        elseif strcmp(usedLabels{1}, 'tantotempo');
            Y.L(labels) = 17;       
        elseif strcmp(usedLabels{1}, 'buonissimo');
            Y.L(labels) = 18;       
        elseif strcmp(usedLabels{1}, 'messidaccordo');
            Y.L(labels) = 19;       
        elseif strcmp(usedLabels{1}, 'sonostufo');
            Y.L(labels) = 20;  
        elseif strcmp(usedLabels{1}, 'idle');
            Y.L(labels) = 21; 
        end            
        labels = strcmp(usedLabels{1}, usedLabels);
        usedLabels(labels) = [];
        i = i + 1;
    end
    if ~all(Y.L)
        error('There is some 0 label. Check the label correspondence to the original data');
    end
else
    if nframesSeg       % random proportional assignment
        if ~nSegs
            error('Number of segments must be greater than 0');
        end
        minCsamples = round(length(toMat(X))*nframesSeg/nSegs);
        n = minCsamples; 
        classes = 1:1:nSegs; L = [];
        c = randi([1 nSegs]); 

        for i = 1:length(X),
            n = n - nframesSeg;
            if n <= 0 || i == length(X),
                n = minCsamples;        
                L = [L c];
                classes(classes==c) = [];
                while ~ismember(c,classes) && ~isempty(classes),
                    c = randi([1 nSegs]); 
                end        
            end
        end
        Y.L = L;
    end
end