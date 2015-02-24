function probs=evaluateSequences(clusters,data,hmmTR, hmmE)
    % Analyze each sequence independently
    if ~iscell(data)
        seq=data;
        n=1;
    else
        n=length(data);
    end
    for i=1:n,
        % Get the sequence
        if iscell(data)
            probs=zeros(1,length(data));
            seq=cell2mat(data(i));
        end
%         if min(seq)<1,
%             display('zero');
%         end
        % Discretize the sequence
        discSeq=discretizeSequence(clusters,seq);
        % Evaluate sequence        
        probs(i) = evaluateHMM(discSeq,hmmTR, hmmE);
    end