function [discreteSeq]=discretizeSequence(c,seq)
    if ~iscell(seq),
        discreteSeq=zeros(size(seq,1),1);
        for i=1:size(seq,1),
            discreteSeq(i)=getCluster(c,seq(i,:));
            if min(discreteSeq(i))<1,
                display(i);
            end
        end
    else
        discreteSeq={};
        for i=1:length(seq),
            discreteSeq=[discreteSeq discretizeSequence(c,seq{i})];
        end
    end
    