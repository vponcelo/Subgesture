function [model,score,bestScores] = evalswHMM(model, X, Y)
% function that returns the score and learn parameters from the evaluation

global DATATYPE; global NAT;

if iscell(model.phmm.model{1}), modelpmtk = model.phmm.model{1}; else modelpmtk = model.phmm.model; end
if iscell(model.phmm.hmmTR_f{1}), TRANS = model.phmm.hmmTR_f{1}; else TRANS = model.phmm.hmmTR_f; end
if iscell(model.phmm.hmmE_f{1}), EMIS = model.phmm.hmmE_f{1}; else EMIS = model.phmm.hmmE_f; end

if iscell(X)
    probs = zeros(1,length(X));
    for i = 1:length(X)
        for j = 1:length(TRANS)
            probs(i,j) = evaluateSequences([],X{i}, TRANS{j}, EMIS{j}, modelpmtk{j});
        end
    end
    %% number of tresholds and interval to test threshs
    
    scoresP = zeros(length(TRANS),model.nThreshs);
    scoresR = zeros(length(TRANS),model.nThreshs);
    scoresA = zeros(length(TRANS),model.nThreshs);
    
    thresholds = zeros(length(TRANS),model.nThreshs);
    
    %% accuracy estimation for single-label prediction
    Yones = ones(1,length(X))';
    for j = 1:length(TRANS)
        interv = (max(probs(:,j))-min(probs(:,j)))/model.nThreshs;
        tMin = min(probs(:,j));
        if interv == 0, interv = tMin*2/model.nThreshs; end
        thresholds(j,:) = tMin + ((1:model.nThreshs)-1)*interv;
        for i = 1:model.nThreshs,
            idxDet = probs(:,j) >= thresholds(j,i);
            TP=sum(idxDet & Yones);            
            TN=sum(~idxDet & ~Yones);
            FP=sum(idxDet & ~Yones);
            FN=sum(~idxDet & Yones);

            if model.accuracyglobal,
                %%% global accuracy:
                scoresA(j,i)= (TN+TP)./(TP+TN+FP+FN);
            else
                %%% weighted accuracy (for imbalanced data sets). 
                scoresA(j,i)= (TP/(TP+FN) + TN/(TN+FP))/2;
            end

            scoresP(j,i) = TP./(TP+FP);
            scoresR(j,i) = TP./(TP+FN);

            if isnan(scoresP(j,i)), scoresP(j,i) = 0; end
            if isnan(scoresR(j,i)), scoresR(j,i) = 0; end
            if isnan(scoresA(j,i)), scoresA(j,i) = 0; end
        end
    end
    %% get the predicted class giving the maximum scores
    k = zeros(1,3); bestScores = zeros(1,3);
    mScoresP = max(scoresP'); mScoresR = max(scoresR'); mScoresA = max(scoresA');
    kscoresP = ismember(mScoresP,max(mScoresP));
    kscoresR = ismember(mScoresR,max(mScoresR)); 
    kscoresA = ismember(mScoresA,max(mScoresA));
    if kscoresP(Y{1}(1)) == 1,
        bestScores(1) = max(mScoresP); k(1) = Y{1}(1);
    else
        [bestScores(1),k(1)] = max(mScoresP);
    end
    if kscoresR(Y{1}(1)) == 1,
        bestScores(2) = max(mScoresR); k(2) = Y{1}(1);
    else
        [bestScores(2),k(2)] = max(mScoresR);
    end
    if kscoresA(Y{1}(1)) == 1,
        bestScores(3) = max(mScoresA); k(3) = Y{1}(1);
    else
        [bestScores(3),k(3)] = max(mScoresA);
    end 
    score = k(model.score2optim);
else    % X evaluate whole seq.
    if model.phmm.pmtk, nm=length(modelpmtk); else nm=length(TRANS); end
    
%     Y.L=Y.Lfr; d=diff(Y.Lfr); Y.L(d==0)=[]; Y.seg=[1 find(d~=0)+1 length(Y.Lfr)];
    predictions = [];
    
    %% number of tresholds and interval to test threshs
    if isempty(model.bestThs), 
        thresholds = zeros(nm,model.nThreshs);
    else
        thresholds(:,1) = model.bestThs; 
        if model.nThreshs > 1, error('g:nThErr','Only 1 threshold should have been learnt per class'); end
    end
    predictions = []; scoresP = zeros(nm,model.nThreshs); scoresR = zeros(nm,model.nThreshs); scoresA = zeros(nm,model.nThreshs);
    
    %% begin evaluation    
    if model.classification
        display(sprintf('Classifying %d sequences in %d gesture classes ...',length(Y.L),nm));
        %% Estimate scores & predictions from probabilities of generating X from each model
        probs = zeros(length(Y.L),nm);
        for s = 1:length(Y.L)      % save each sequence vs model dtw costs
            if s < length(Y.L)
                seq = X(Y.seg(s):Y.seg(s+1)-1);
            else
                seq = X(Y.seg(s):Y.seg(s+1));
            end
            for k = 1:length(TRANS)
                probs(s,k) = evaluateSequences([], seq, TRANS{k}, EMIS{k}, modelpmtk{k});
            end
        end
        if isempty(model.bestThs)
            for k = 1:nm     % save dtw costs of each non-iddle sequence vs its model
                interv = (max(probs(Y.L==k,k))-min(probs(Y.L==k,k)))/model.nThreshs;
                tMin = min(probs(Y.L==k,k));
                if interv == 0, interv = tMin*2/model.nThreshs; end
                thresholds(k,:) = tMin + ((1:model.nThreshs)-1)*interv;
            end
        end
        %% accuracy estimation for single-label prediction 
        Yones=zeros(size(Y.L,2),nm);
        for i=1:nm,
            Yones(Y.L==i,i)=1;
        end
        for i = 1:model.nThreshs,
            for s = 1:nm,
                idxDet = probs(:,s) >= thresholds(s,i);  
                TP=sum(idxDet & Yones(:,s));            
                TN=sum(~idxDet & ~Yones(:,s));
                FP=sum(idxDet & ~Yones(:,s));
                FN=sum(~idxDet & Yones(:,s));

                if model.accuracyglobal,
                    %%% global accuracy:
                    scoresA(s,i)= (TN+TP)./(TP+TN+FP+FN);
                else
                    %%% weighted accuracy (for imbalanced data sets). 
                    scoresA(s,i)= (TP/(TP+FN) + TN/(TN+FP))/2;
                end

                scoresP(s,i) = TP./(TP+FP);
                scoresR(s,i) = TP./(TP+FN);

                if isnan(scoresP(s,i)), scoresP(s,i) = 0; end
                if isnan(scoresR(s,i)), scoresR(s,i) = 0; end
                if isnan(scoresA(s,i)), scoresA(s,i) = 0; end
            end
        end
        [bestScores(1),bestThsPos(1)] = max(mean(scoresP));
        [bestScores(2),bestThsPos(2)] = max(mean(scoresR));
        [bestScores(3),bestThsPos(3)] = max(mean(scoresA));
        score = bestScores(model.score2optim);
        model.bestThs = zeros(1,nm); model.nThreshs = 1;
        for i = 1:nm
            model.bestThs(i) = thresholds(i,bestThsPos(model.score2optim));
        end
    else
        %% Spotting
        display(sprintf('Spotting the %d classes in a sequence of %d frames ...',nm,length(Y.Lfr)));
        global DATATYPE; global NAT;
        sws = 10:10:60;     % size of gesture sequences is usually < 60
        
        %% number of tresholds and interval to test threshs
        if ~isempty(model.bestThs)
            if model.nThreshs > 1, error('g:nThErr','Only 1 threshold should have been learnt per class'); end
            thresholds(:,1) = model.bestThs;
        else
            thresholds = zeros(nm,model.nThreshs);
        end
        
        scoresO = zeros(nm,model.nThreshs); scoresP = zeros(nm,model.nThreshs); scoresR = zeros(nm,model.nThreshs); scoresA = zeros(nm,model.nThreshs);
        prexp = false(model.nThreshs,length(Y.Lfr),nm);
        bestScores = zeros(nm,4); bestThsPos = zeros(nm,4);
        %% Estimate scores & predictions from probabilities of generating X from each model
        for k = 1:nm
            GTtestk = Y.L == k; GTtestkFr = Y.Lfr == k;
            probs = cell(1,length(sws));
            for i = 1:length(sws)
                probs{i} = [];
                for st = 1:sws(i):size(X,2)-1
                    seq = X(:,st:min(st+sws(i)-1,size(X,2)));
                    probs{i} = [probs{i} evaluateSequences([], seq, TRANS{k}, EMIS{k}, modelpmtk{k})];
                end
            end
            if isempty(model.bestThs)
                maxim = max(cellfun(@max,probs)); minim = min(cellfun(@min,probs));
                interv = (maxim-minim)/model.nThreshs; tMin = minim;
                if interv == 0, interv = tMin*2/model.nThreshs; end
                thresholds(k,:) = tMin + ((1:model.nThreshs)-1)*interv;
            end
            for j = 1:model.nThreshs
                for i = 1:length(sws)
                    det = probs{i} >= thresholds(k,j);
                    det = reshape(repmat(det,sws(i),1),[1 size(det,2)*sws(i)]);  % replicate to obtain detected frames
                    if length(det) > size(X,2)
                        det(size(X,2)+1:end) = [];                               % remove offset sw
                    else
                        det = [det repmat(det(end),size(X,2)-length(det))];      % add offset sw
                    end
                    prexp(j,:,k) = prexp(j,:,k) | det;
                end
                if strcmp(DATATYPE,'chalearn2014') && NAT == 3
                    prexp(j,:,k) =([prexp(j,6:end,k),0,0,0,0,0]);       % correct offset of deep features
                end                
                %% compute scores from evaluated probabilities and thresholds
                detSw = prexp(j,:,k);
                scoresO(k,j) = sum(GTtestkFr & detSw)./sum(GTtestkFr | detSw);     % overlap (Jaccard Index)

                detSw = getActivations(detSw, GTtestkFr, Y.seg, model);
                % only for MADX database (recognition)
                if strcmp(DATATYPE,'mad1') || strcmp(DATATYPE,'mad2') ...
                        || strcmp(DATATYPE,'mad3') || strcmp(DATATYPE,'mad4') ...
                        || strcmp(DATATYPE,'mad5') 
                    [~,~,R] = estimate_overlap_madold(GTtestk, detSw, model.minOverlap);
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
            [bestScores(k,1),bestThsPos(k,1)] = max(scoresO(k,:));
            [bestScores(k,2),bestThsPos(k,2)] = max(scoresP(k,:));
            [bestScores(k,3),bestThsPos(k,3)] = max(scoresR(k,:));
            [bestScores(k,4),bestThsPos(k,4)] = max(scoresA(k,:));            
            %% save mean scores and learnt thresholds (allow multi-label assignment per frame)  
%            score = mean(bestScores(:,model.score2optim));
%             bestScores = mean(bestScores)
            
%            gtF = Y.Lfr;
%            save('predictions.mat','predictionsF','gtF');
%            imagesc(confusionmat(predictionsF,gtF)); colormap(hot);
        
            if strcmp(model.scoreMeasure,'levenshtein')
                [~,pos] = max(scoresO{k});
                idx = find(prexp(pos,:,k) == 1);
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
        
        if strcmp(model.scoreMeasure,'overlap')            
            %% here we generate a single prediction per frame
            predictions=-ones(1,size(prexp,2)); 
            for i=1:size(prexp,2),       
               cupred= zeros(1,k);
               cues = reshape(prexp(:,i,:),size(prexp,1),size(prexp,3));
               for j=1:k,
                   cupred(j)=cues(bestThsPos(j,1),j);           
               end

               ofin=find(cupred~=0);
               if ~isempty(ofin),           
                   if length(ofin)==1
                       predictions(i)=ofin; % a single gesture activated
                   else
        %                predictions(i)=ofin(randi(length(ofin))); % random prediction               
                        ofinscores=bestScores(ofin,model.score2optim); % based on performance
                        [~,sb]=max(ofinscores); predictions(i)=ofin(sb);
                   end
               end       
            end
%             mean(bestScores)
            clear bestScores;

            %%% if we want to plot predictions vs GT
            % plot(predictions);
            % gcf;hold on;
            % plot(nY,':r');

            % now estimate overlap, prec, rec, and f1 (no accuracy right now)
            [~,~,R,ovlp] = estimate_overlap_mad(Y.Lfr,predictions,model.minOverlap);

            % Assign Overlap, Precision, Recall, F1-score
            bestScores(1) = ovlp;  bestScores(2) = R.prec; bestScores(3) = R.rec; bestScores(4) = (2.*R.rec.*R.prec)./(R.rec + R.prec);
%             bestScores
            for j=1:length(bestScores), if isnan(bestScores(j)), bestScores(j)=0; end; end
            
            switch model.score2optim
                case 1, score=ovlp;
                case 2, score=R.prec;
                case 3, score=R.rec;
                case 4, score=(2.*R.rec.*R.prec)./(R.rec + R.prec);
            end
            clear prexp R ovlp ofinscores ofin cupred cues;
            if isempty(model.bestThs)
                model.nThreshs = 1; model.bestThs = zeros(1,k);
                for i = 1:k, model.bestThs(i) = thresholds(i,bestThsPos(i,model.score2optim)); end
            end
        elseif strcmp(model.scoreMeasure,'levenshtein') && ~isempty(predictions)
            if ~isempty(predictions), plotmistakes(predictions,Y,1); end
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