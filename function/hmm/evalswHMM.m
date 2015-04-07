function [model,score] = evalswHMM(model, X, Y, TRANS, EMIS, modelpmtk)
% function that returns the score and learn parameters from the evaluation

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
    if model.classification
        Y.L=Y.Lfr; d=diff(Y.Lfr); Y.L(d==0)=[]; Y.seg=[1 find(d~=0)+1];
    end
    if model.phmm.pmtk, nm=length(modelpmtk); else nm=length(TRANS); end
    swScores = cell(1,length(nm)); 
    thresholds = cell(1,length(nm));
    % test each gesture
    for k = 1:nm
        if model.classification
            GTtestk = Y.L == k;            
        else
            GTtestk = Y.Lfr == k;
        end
        probs = cell(1,length(sws));
        %% evaluate sliding windows
        for i = 1:length(sws)
            probs{i} = [];
            for st = 1:sws(i):size(X,2)-1
                seq = X(:,st:min(st+sws(i),size(X,2)));
                probs{i} = [probs{i} evaluateSequences([], seq, TRANS{k}, EMIS{k}, modelpmtk{k})];
            end
        end
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
        swScores{k} = zeros(length(sws),model.nThreshs);
        %% compute scores from evaluated probabilities
        thresholds{k} = zeros(1,model.nThreshs);
        for i = 1:length(sws)
            for j = 1:model.nThreshs
                if isempty(model.bestThs)
                    thresholds{k}(j) = tMin + (j-1)*interv;
                else
                    thresholds{k}(j) = model.bestThs(k);
                end
                detSw = probs{i} > thresholds{k}(j); detSw = reshape(repmat(detSw,sws(i),1),[1 size(detSw,2)*sws(i)]);
                if length(detSw) > size(X,2)
                    detSw(size(X,2)+1:end) = [];                                    % remove offset sw
                else
                    detSw = [detSw repmat(detSw(end),size(X,2)-length(detSw))];     % add offset sw
                end
                detSw=([detSw(6:end),0,0,0,0,0]);
                if model.classification
                    d = diff(detSw); idxL = find(d~=0); 
                    detSw = zeros(1,length(GTtestk));
                    for l = 1:length(idxL)
                        cLabel = sum(Y.seg < idxL(l));
                        if cLabel > 0
                            detSw(cLabel) = 1;
                        end
                    end
                end
                if model.classification
                    swScores{k}(i,j) = sum(GTtestk & detSw)./sum(GTtestk & detSw | ~GTtestk & detSw);   % Precision
                    swScores{k}(i,j) = sum(GTtestk & detSw)./sum(GTtestk & detSw | GTtestk & ~detSw);  % Recall
                    swScores{k}(i,j) = sum(GTtestk & detSw | ~GTtestk & ~detSw)./sum(GTtestk & detSw | GTtestk & ~detSw | ~GTtestk & detSw | ~GTtestk & ~detSw);   % Accuracy
                else
                    swScores{k}(i,j) = sum(GTtestk & detSw)./sum(GTtestk | detSw);          % overlap (Jaccard Index)
                end
            end
        end
    end
    bestScores = cell(1,length(sws)); bestThsPos = cell(1,length(sws));
    score = zeros(1,length(sws));
    for i = 1:length(sws)
        bestScores{i} = zeros(1,k); bestThsPos{i} = zeros(1,k);
        for l = 1:k
            [bestScores{i}(l),bestThsPos{i}(l)] = max(swScores{k}(i,:));
        end
        score(i) = mean(bestScores{i});
    end
    [score,p] = max(score);    
    for i = 1:k
        model.bestThs(i) = thresholds{i}(bestThsPos{p}(i));
    end
end