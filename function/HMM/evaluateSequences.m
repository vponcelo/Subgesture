function probs=evaluateSequences(clusters,data,hmmTR, hmmE)
    % Analyze each sequence independently
    probs=zeros(1,length(data));
    for i=1:length(data),
        % Get the sequence
        seq=cell2mat(data(i));
        if min(seq)<1,
            display('zero');
        end
        % Discretize the sequence
        discSeq=discretizeSequence(clusters,seq);
        % Evaluate sequence        
        probs(i) = evaluateHMM(discSeq,hmmTR, hmmE);
    end