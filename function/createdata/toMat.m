function [D]=toMat(seqs)
    if iscell(seqs)
        D=[];
        for s=seqs,
            D=[D;cell2mat(s)];
        end
    else
        D = seqs;
    end