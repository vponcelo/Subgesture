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
else
    % X evaluate whole seq.
    if model.phmm.pmtk, nm=length(modelpmtk); else nm=length(TRANS); end
    if model.classification
        Y.L=Y.Lfr; d=diff(Y.Lfr); Y.L(d==0)=[]; Y.seg=[1 find(d~=0)+1 length(Y.Lfr)];
        scoresP = cell(1,nm); scoresR = cell(1,nm); scoresA = cell(1,nm);
    end    
    scoresO = cell(1,nm); 
    thresholds = cell(1,length(nm));
    display(sprintf('Evaluating sequence of %d frames ...',length(Y.Lfr)));
    % test each gesture
    for k = 1:nm
        GTtestk = Y.L == k;
        GTtestkFr = Y.Lfr == k;
        %% number of tresholds and interval to test threshs
        if ~isempty(model.bestThs)
            model.nThreshs = 1;
        else
            TOL_THRESH = 0.01; tMin = 0;
            interv = 1./model.nThreshs;
            while interv < TOL_THRESH && model.nThreshs > 1
                model.nThreshs = round(model.nThreshs/2);
                interv = 1./model.nThreshs;
                display(sprintf('Decreasing number of test thresholds to %d',model.nThreshs));
            end
        end
        %% compute scores from evaluated probabilities and thresholds
        thresholds{k} = zeros(1,model.nThreshs);
        detSeqLog = zeros(model.nThreshs,length(GTtestkFr));
        scoresO{k} = zeros(1,model.nThreshs); 
        if model.classification
            scoresP{k} = zeros(1,model.nThreshs); scoresR{k} = zeros(1,model.nThreshs); scoresA{k} = zeros(1,model.nThreshs);
        end
        for j = 1:model.nThreshs
            if isempty(model.bestThs)
                thresholds{k}(j) = tMin + (j-1)*interv;
            else
                thresholds{k}(j) = model.bestThs(k);
            end  
            for i = 1:length(sws)
                probs = [];
                for st = 1:sws(i):size(X,2)-1
                    seq = X(:,st:min(st+sws(i),size(X,2)));
                    probs = [probs evaluateSequences([], seq, TRANS{k}, EMIS{k}, modelpmtk{k})];
                end
                det = probs > thresholds{k}(j);
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
            scoresO{k}(j) = sum(GTtestkFr & detSeqLog(j,:))./sum(GTtestkFr | detSeqLog(j,:));     % overlap (Jaccard Index)
            if isnan(scoresO{k}(j)), scoresO{k}(j) = 0; end
            if model.classification
                % recognition from spotting
                detSw = getActivations(detSeqLog(j,:), GTtestkFr, Y.seg, model);
                % only for MADX database (recognition)
                if strcmp(DATATYPE,'mad1') || strcmp(DATATYPE,'mad2') ...
                        || strcmp(DATATYPE,'mad3') || strcmp(DATATYPE,'mad4') ...
                        || strcmp(DATATYPE,'mad5') 
                    [~,~,R] = estimate_overlap_mad(GTtestk, detSeqLog(j,:), model.minOverlap);
                    scoresP{k}(j) = R.prec2;    % Precision
                    scoresR{k}(j) = R.rec2;     % Recall
                else
                    scoresP{k}(j) = sum(GTtestk & detSw)./sum(GTtestk & detSw | ~GTtestk & detSw);  % Precision
                    scoresR{k}(j) = sum(GTtestk & detSw)./sum(GTtestk & detSw | GTtestk & ~detSw);  % Recall
                end
                scoresA{k}(j) = sum(GTtestk & detSw | ~GTtestk & ~detSw)./sum(GTtestk & detSw | GTtestk & ~detSw | ~GTtestk & detSw | ~GTtestk & ~detSw);   % Accuracy
                if isnan(scoresP{k}(j)), scoresP{k}(j) = 0; end
                if isnan(scoresR{k}(j)), scoresR{k}(j) = 0; end
                if isnan(scoresA{k}(j)), scoresA{k}(j) = 0; end
            end            
        end
    end
    if ~model.classification,
        bestScores = zeros(k,1); bestThsPos = zeros(k,1);
    else
        bestScores = zeros(k,4); bestThsPos = zeros(k,4);
        
    end
    model.bestThs = zeros(1,k);
    for i = 1:k
        [bestScores(i,1),bestThsPos(i,1)] = max(scoresO{i});
        if model.classification
            [bestScores(i,2),bestThsPos(i,2)] = max(scoresP{i});
            [bestScores(i,3),bestThsPos(i,3)] = max(scoresR{i});
            [bestScores(i,4),bestThsPos(i,4)] = max(scoresA{i});            
        end
        if strcmp(model.scoreMeasure,'overlap')
            model.bestThs(i) = max(thresholds{i}(bestThsPos(i,model.score2Optim)));
        end
    end
    end
    if ~model.classification
        score = mean(bestScores(:,1));
    else
        score = mean(bestScores(:,model.score2Optim));     
    end
    bestScores = reshape(mean(bestScores),[1 size(bestScores,2)]);
end