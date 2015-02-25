function probs=evaluateSequences(clusters,data,hmmTR, hmmE)
    % Analyze each sequence independently
    if ~iscell(data)
        n=1;
        seq=data;
    else
        n=length(data);
    end
    probs=zeros(1,n);
    for i=1:n,
        % Get the ith sequence
        if iscell(data)
            seq=cell2mat(data(i));
        end
        if min(seq)<1,
            warning('evaluateSequences:zero','there is a zero in the sequence');
        end
        % Discretize the sequence
        if ~isempty(clusters)
            discSeq=discretizeSequence(clusters,seq);
        else
            discSeq=seq;
        end
        % Evaluate sequence        
        probs(i) = evaluateHMM(discSeq,hmmTR, hmmE);
    end