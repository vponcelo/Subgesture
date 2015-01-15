function [Xseq,Y] = setLabels(X,skels,nframesSeg)
% Set the labels of the sequence either with the ground truth or with a
% fixed number of frames for each segment.
%
% output
%     Xseq: Data sequences
%     Y: Data segments and labels
% input
%     X: whole data in frame list
%     skels: skeleton sequences
%     nframesSeg: fixed number of fixed frames to segment the sequences

if ~nframesSeg
    Y.cnames = []; Y.seg = 1; iseg = 1; totalSkels = 0;
    for i = 1:length(skels)        
        totalSkels = totalSkels + length(skels{i}.skeleton);
        for j = 1:length(skels{i}.labels)
            while iseg < skels{i}.labels(j).End
                if skels{i}.labels(j).Begin-1 > iseg
                    Y.seg = [Y.seg Y.seg(end)-iseg+skels{i}.labels(j).Begin]; 
                    Y.cnames = [Y.cnames {'idle'}];
                    iseg = skels{i}.labels(j).Begin;
                elseif iseg == skels{i}.labels(j).Begin
                    Y.seg = [Y.seg Y.seg(end)+skels{i}.labels(j).End-iseg+1];
                    Y.cnames = [Y.cnames {skels{i}.labels(j).Name}];
                    iseg = skels{i}.labels(j).End+1;
                else
                    iseg = iseg + 1;
                end
            end
            if j == length(skels{i}.labels) 
                if Y.seg(end) == totalSkels          
                    Y.seg(end) = totalSkels+1;                    
                elseif Y.seg(end) < totalSkels
                    Y.seg = [Y.seg totalSkels+1];
                    Y.cnames = [Y.cnames {'idle'}];
                end
                if i == length(skels) && Y.seg(end) > length(X)
                    Y.seg(end) = Y.seg(end)-1;
                end
            end
        end
        iseg = 1;
    end
    Xseq = X;
else
    Y = [];
    nframes = size(X,1);
    Xseq = cell(1,800);
    i = 1; in = i;
    while in < nframes,
        en = in + nframesSeg - 1;
        if en >= nframes, 
            en = nframes;
        end
        Xseq{i} = X(in:en,:);  
        in = en + 1;
        i = i + 1;
    end
    Xseq(i:end) = [];

    for i = 1:length(Xseq),
        M = [];
        for j = 1:size(Xseq,1),
            M = horzcat(M,Xseq{i}(j,:));
        end
        Xseq{i} = M;
    end
end

