function s = evalswHMM(seqTest, TRANS, EMIS, modelpmtk)
% function that returns the maximum probability of the evaluation

sws = 10:10:60;     % size of gesture sequences is usually < 60
if iscell(seqTest)
    scores = zeros(length(sws),length(seqTest));
    seqs = cell(1,length(seqTest));
    for i = 1:length(sws)    
        for j = 1:length(seqs)
            seqs{j} = seqTest{j}(:,1:min(sws(i),size(seqTest{j},2)));
        end
        scores(i,:) = evaluateSequences([],seqs, TRANS, EMIS, modelpmtk);
    end
    s = max(scores);
else
    % seqTest evaluate whole seq.
    s = []; st = 1;
    while st + sws(end) <= length(seqTest)
        scores = zeros(length(sws),length(seqTest));
        for i = 1:length(sws)    
            seqs = seqTest(:,st:st+min(sws(i),size(seqTest,2)));
            scores(i,:) = evaluateSequences([],seqs, TRANS, EMIS, modelpmtk);
        end
        s = [s max(scores)];
        st = st + sws(1) + 1;
    end
    s = mean(s);
end