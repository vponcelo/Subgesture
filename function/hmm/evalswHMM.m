function [model,score,bestScores] = evalswHMM(model, X, Y, TRANS, EMIS, modelpmtk)
% function that returns the score and learn parameters from the evaluation

global DATATYPE; global NAT;

sws = 10:10:60;     % size of gesture sequences is usually < 60
if iscell(X)
    probs = zeros(length(sws),length(X));
    seqs = cell(1,length(X));
    for i = 1:length(sws)
        for j = 1:length(seqs)
            seqs{j} = X{j}(:,1:min(sws(i),size(X{j},2)));
        end
        probs(i,:) = evaluateSequences([],seqs, TRANS, EMIS, modelpmtk);
    end
    score = max(probs); model = [];
else    % X evaluate whole seq.
    if model.phmm.pmtk, nm=length(modelpmtk); else nm=length(TRANS); end
    
    Y.L=Y.Lfr; d=diff(Y.Lfr); Y.L(d==0)=[]; Y.seg=[1 find(d~=0)+1 length(Y.Lfr)];
    predictions = [];
    
    %% number of tresholds and interval to test threshs
    if ~isempty(model.bestThs)
        model.nThreshs = 1;
        thresholds(:,1) = model.bestThs;
    else
        thresholds = zeros(nm,model.nThreshs);
        interv = 1./model.nThreshs;
        for k = 1:nm     % save dtw costs of each non-iddle sequence vs its model
            thresholds(k,:) = (1:model.nThreshs)*interv;
        end
    end
    
    %% begin evaluation
    display(sprintf('Evaluating sequence of %d frames in %d gesture classes ...',length(Y.Lfr),nm));    
    if model.classification        
        probs = zeros(length(Y.L),nm);
        for s = 1:length(Y.L)      % save each sequence vs model dtw costs
            if s < length(Y.L)
                seq = X(Y.seg(s):Y.seg(s+1)-1);
            else
                seq = X(Y.seg(s):Y.seg(s+1));
            end
            for k = 1:nm
                probs(s,k) = evaluateSequences([], seq, TRANS{k}, EMIS{k}, modelpmtk{k});                
            end
        end
        
        scoresP = zeros(length(Y.L),model.nThreshs); scoresR = zeros(length(Y.L),model.nThreshs); scoresA = zeros(length(Y.L),model.nThreshs);
        for s = 1:length(Y.L)
            for i = 1:model.nThreshs
                % Se usan las 
                TP = 0; FP = 0; FN = 0; %TN = 0;    % NO SE CONSIDERAN LOS TN, DECIDIR SI LOS USAMOS SEGUN LO HABLADO
                idxDet = probs(s,:) < thresholds(:,i)';
                if Y.L(s) < nm+1    % iddle gesture is ignored if it wasn't learnt
                    if idxDet(Y.L(s))
                        TP = TP + 1;    % DECIDISION DE LA ASIGNACION DE VERDADEROS POSITIVOS
                        if sum(idxDet) > 1, FP = FP + sum(idxDet)-1; end
                    else
                        FN = FN + 1; FP = FP + sum(idxDet); 
                    end
                else
                    FP = FP + sum(idxDet);
                end
%                 if sum(~idxDet)
%                     if ~idxDet(Y.L(s)) 
%                         TN = TN + sum(~idxDet)-1;
%                     else
%                         TN = TN + sum(~idxDet); 
%                     end
%                 end
                scoresP(s,i) = TP./(TP+FP);
                scoresR(s,i) = TP./(TP+FN);
                scoresA(s,i) = (TP)./(TP+FN+FP);    %%% (TP/(TP+FN) + TN/(TN+FP))/2   % Métrica para el balanceo entre clases
                if isnan(scoresP(s,i)), scoresP(s,i) = 0; end
                if isnan(scoresR(s,i)), scoresR(s,i) = 0; end
                if isnan(scoresA(s,i)), scoresA(s,i) = 0; end
            end
        end
        [bestScores(1),bestThsPos(1)] = max(mean(scoresP));
        [bestScores(2),bestThsPos(2)] = max(mean(scoresR));
        [bestScores(3),bestThsPos(3)] = max(mean(scoresA));
        score = bestScores(model.score2optim);
        model.bestThs = zeros(1,k);
        for i = 1:k
            model.bestThs(i) = thresholds(i,bestThsPos(model.score2optim));
        end
    else
        scoresO = zeros(nm,model.nThreshs); scoresP = zeros(nm,model.nThreshs); scoresR = zeros(nm,model.nThreshs); scoresA = zeros(nm,model.nThreshs);
        thresholds = zeros(nm,model.nThreshs);        
            
        display(sprintf('Evaluating sequence of %d frames in %d gesture classes ...',length(Y.Lfr),nm));
        % test each gesture
        for k = 1:nm
            GTtestk = Y.L == k;
            GTtestkFr = Y.Lfr == k;
            %% compute scores from evaluated probabilities and thresholds
            detSeqLog = zeros(model.nThreshs,length(GTtestkFr));
            for j = 1:model.nThreshs
                if isempty(model.bestThs)
                    thresholds(k,j) = tMin + (j-1)*interv;
                else
                    thresholds(k,j) = model.bestThs(k);
                end
                for i = 1:length(sws)
                    probs = [];
                    for st = 1:sws(i):size(X,2)-1
                        seq = X(:,st:min(st+sws(i),size(X,2)));
                        probs = [probs evaluateSequences([], seq, TRANS{k}, EMIS{k}, modelpmtk{k})];
                    end
                    det = probs > thresholds(k,j);
                    det = reshape(repmat(det,sws(i),1),[1 size(det,2)*sws(i)]);  % replicate to obtain detected frames
                    if length(det) > size(X,2)
                        det(size(X,2)+1:end) = [];                               % remove offset sw
                    else
                        det = [det repmat(det(end),size(X,2)-length(det))];      % add offset sw
                    end
                    detSeqLog(j,:) = detSeqLog(j,:) | det;
                end
                if strcmp(DATATYPE,'chalearn2014') && NAT == 3
                    detSeqLog(j,:) =([detSeqLog(j,6:end),0,0,0,0,0]);       % correct offset of deep features
                end
                scoresO(k,j) = sum(GTtestkFr & detSeqLog(j,:))./sum(GTtestkFr | detSeqLog(j,:));     % overlap (Jaccard Index)
                if model.classification
                    % recognition from spotting
                    detSw = getActivations(detSeqLog(j,:), GTtestkFr, Y.seg, model);
                    % only for MADX database (recognition)
                    if strcmp(DATATYPE,'mad1') || strcmp(DATATYPE,'mad2') ...
                            || strcmp(DATATYPE,'mad3') || strcmp(DATATYPE,'mad4') ...
                            || strcmp(DATATYPE,'mad5') 
                        [~,~,R] = estimate_overlap_mad(GTtestk, detSeqLog(j,:), model.minOverlap);
                        scoresP(k,j) = R.prec2;    % Precision
                        scoresR(k,j) = R.rec2;     % Recall
                    else
                        scoresP(k,j) = sum(GTtestk & detSw)./sum(GTtestk & detSw | ~GTtestk & detSw);  % Precision
                        scoresR(k,j) = sum(GTtestk & detSw)./sum(GTtestk & detSw | GTtestk & ~detSw);  % Recall
                    end
                    scoresA(k,j) = sum(GTtestk & detSw | ~GTtestk & ~detSw)./sum(GTtestk & detSw | GTtestk & ~detSw | ~GTtestk & detSw | ~GTtestk & ~detSw);   % Accuracy
                    if isnan(scoresO(k,j)), scoresO(k,j) = 0; end
                    if isnan(scoresP(k,j)), scoresP(k,j) = 0; end
                    if isnan(scoresR(k,j)), scoresR(k,j) = 0; end
                    if isnan(scoresA(k,j)), scoresA(k,j) = 0; end
                end
            end
            if strcmp(model.scoreMeasure,'levenshtein')
                [~,pos] = max(scoresO(k,:));
                idx = find(detSeqLog(pos,:) == 1);
                if ~isempty(idx)
                    inF = idx(1);
                    for i = 2:length(idx)
                        if idx(i)-idx(i-1) > 1
                            endF = idx(i-1);
                            predictions = [predictions inF endF k];
                            inF = idx(i);
                        end
                        if i == length(idx)
                            endF = idx(i);
                            predictions = [predictions inF endF k];
                            % Last ones when inF == endF will be also added.
                        end
                    end
                else
                    predictions = 0;
                end
                predLabels = []; lints = [];
            end
        end
        bestScores = zeros(k,4); bestThsPos = zeros(k,4);
        model.bestThs = zeros(1,k);
        for i = 1:k
            [bestScores(i,1),bestThsPos(i,1)] = max(scoresO(i,:));
            [bestScores(i,2),bestThsPos(i,2)] = max(scoresP(i,:));
            [bestScores(i,3),bestThsPos(i,3)] = max(scoresR(i,:));
            [bestScores(i,4),bestThsPos(i,4)] = max(scoresA(i,:));
            model.bestThs = thresholds(i,bestThsPos(i,model.score2optim));
        end
        score = mean(bestScores(:,model.score2optim));
        bestScores = mean(bestScores);
        
        if ~isempty(predictions), plotmistakes(predictions,Y,1); display(e.message); end
        if strcmp(model.scoreMeasure,'levenshtein') && ~isempty(predictions)
            for i = 1:size(X,1)
                l = find(predictions(1:3:end) == i);
                while ~isempty(l) && ~any(ismember(l,lints))
                    [minFi,pos] = min(predictions(3*(l-1)+2));
                    predLabels = [predLabels predictions(3*l(pos))];
                    l(pos) = [];
                    lint = find(predictions(1:3:end) > i & predictions(2:3:end) < minFi);
                    lint(ismember(lint,lints)) = [];
                    lints = [lints lint];
                    while ~isempty(lint)
                        [~,posInt] = min(predictions(3*(lint-1)+1));
                        predLabels = [predLabels predictions(3*lint(posInt))];
                        lint(posInt) = [];
                    end
                end
            end
            score = levenshtein(predLabels,Y.L);
        end
    end
end