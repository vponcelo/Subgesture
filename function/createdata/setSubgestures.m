function [Xseq,Y] = setSubgestures(X,skels,nframesSeq)
% output
%     Xseq: Data sequences
%     Y: Data segments and labels
% input
%     X: whole data in frame list
%     skels: skeleton sequences
%     nframesSeq: fixed number of fixed frames per subgesture sequence (set for 1 non subgesture)

if ~nframesSeq
    Y.cnames = []; Y.seg = 1; k = 0; iseg = 1;
    for i = 1:length(skels)        
        for j = 1:length(skels{i}.labels)
            while iseg < skels{i}.labels(j).End
                if skels{i}.labels(j).Begin-1 > iseg
                    Y.seg = [Y.seg Y.seg(end)-iseg+skels{i}.labels(j).Begin-1]; 
                    Y.cnames = [Y.cnames {'idle'}];
                    iseg = skels{i}.labels(j).Begin;
                elseif iseg == skels{i}.labels(j).Begin
                    Y.seg = [Y.seg Y.seg(end)+skels{i}.labels(j).End-iseg+1];
                    Y.cnames = [Y.cnames {skels{i}.labels(j).Name}];
                    iseg = skels{i}.labels(j).End;
                else
                    iseg = iseg + 1;
                end
            end
            if j == length(skels{i}.labels) 
                if ~k 
                    Y.seg = [Y.seg length(skels{i}.skeleton)];                     
                else
                    Y.seg = [Y.seg Y.seg(k)+length(skels{i}.skeleton)];                      
                end
                k = length(Y.seg);
                Y.cnames = [Y.cnames {'idle'}];
            end
        end
        iseg = 0;
    end
    Xseq = X;
else
    Y = [];
    nframes = size(X,1);
    Xseq = cell(1,800);
    i = 1; in = i;
    while in < nframes,
        en = in + nframesSeq - 1;
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

