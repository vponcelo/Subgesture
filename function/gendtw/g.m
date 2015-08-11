function [model,score,bestScores,predictions] = g(model,Xc,Yc)
% Output:
%   model
%   score(default): mean overlap for the k detected gestures
% Input:
%   model: model parameters
%   Xc: data to test against (cell or whole sequence)
%   Yc: labels of X (cell or whole sequence)

if iscell(Xc) && ~model.classification
    error('g:classErr','Input data cannot be a cell here, need to specify this functionality');
elseif ~model.classification
    sw = model.sw;
    if sw > length(Yc.Lfr)
        error('g:swErr','sliding window is larger than the validation sequence');
    end

    if sw == 0
        sw = length(Xc)-1;
        r = 1;
    else
        sw = model.sw;
        r = inf;
        while any(r > length(Yc.Lfr)-sw)
            r = randperm(round(length(Yc.Lfr)),1);
        end
    end
end

%% select sequence subset to evaluate
nm = length(model.M);
if model.sw == 0,
    X = Xc; Y = Yc;
else
    if ~iscell(Xc)
        seg=r:min(r+sw,length(Yc.Lfr));
        X=Xc(seg,:); Y.Lfr=Yc.Lfr(seg);
        idxSeg = find(Yc.seg >= seg(1) & Yc.seg <= seg(end));
        Y.seg = Yc.seg(idxSeg); Y.L = Yc.L(idxSeg(1:end-1));
        if Y.seg(1) > seg(1), Y.seg = [seg(1) Y.seg]; Y.L = [Yc.L(idxSeg(1)-1) Y.L]; end
        if Y.seg(end) < seg(end), Y.seg = [Y.seg seg(end)]; Y.L = [Y.L Yc.L(idxSeg(end))]; end
        % IMPLEMENT: generate Y.seg form Y.L. This is not necessary for now,
        % just in case we might need a sliding window of size model.sw > 0
    end
end

 %% number of tresholds and interval to test threshs
if isempty(model.bestThs), 
    thresholds = zeros(nm,model.nThreshs);
else
    thresholds(:,1) = model.bestThs; 
    if model.nThreshs > 1, error('g:nThErr','Only 1 threshold should have been learnt per class'); end
end
predictions = []; scoresP = zeros(nm,model.nThreshs); scoresR = zeros(nm,model.nThreshs); scoresA = zeros(nm,model.nThreshs);
    
% begin evaluation
%% classification
if model.classification
    if iscell(X)
        %% begin evaluation
        display(sprintf('Classifying sequences in %d gesture classes ...',nm));
        Wc = cell(1,nm);
        for l = 1:nm
            Wc{l} = zeros(length(X{l}),nm);
            for s = 1:length(X{l})
                seq = X{l}{s};
                for k = 1:nm
                    M = model.M{k};
                    if ~isempty(model.D)
                        if model.darwin
                            KT = model.KT{l}{s};
                        else
                            KT = getUpdatedCosts(seq,model.SM);
                        end
                        if ~iscell(M)
                            W = single(dtwc(seq,M,true,Inf,model.D,model.KM{k},KT));
                        else
                            if ~isempty(M{model.k})
                                W = single(dtwc(seq,M{model.k},true,Inf,model.D,model.KM{k},KT));
                            end
                        end
                    else
                        if ~iscell(M)
                            if ~model.pdtw,
                                W = single(dtwc(seq,M,true));
                            else
                                Pql=zeros(size(seq,1),size(M,1));
                                for hh=1:size(M,1),
                                    if ~isempty(model.lmodel(k,hh).obj),
                                        Pql(:,hh)= mixGaussLogprob(model.lmodel(k,hh).obj,seq);
            %                             Pql(:,hh)= log(pdf(model.lmodel(k,hh).obj,seq));
                                        Pql(isinf(Pql(:,hh)))=0;
                                    end
                                end
                                noze = sum(Pql)~=0;
                                Pql(Pql==0) = mean(mean(Pql(:,noze)));
                                maval = max(max(abs(Pql)));
            %                     Dima  = (1-Pql)./maval;
            %                     DD=pdist2(seq,M);
            %                     DD = DD./max(max(DD));
                                Dima = (1-Pql)./maval.*pdist2(seq,M);
                                W = single(dtw3(seq,M,true,Inf,Dima));
                            end
                        else
                            if ~isempty(M{model.k})
                                W = single(dtwc(seq,M{model.k},true));
                            end
                        end
                    end
                    Wc{l}(s,k) = W(end,end);
                end
            end
            if isempty(model.bestThs)
                interv = (max(Wc{l}(:,l))-min(Wc{l}(:,l)))/model.nThreshs;
                tMin = min(Wc{l}(:,l));
                if interv == 0, interv = tMin*2/model.nThreshs; end
                thresholds(l,:) = tMin + ((1:model.nThreshs)-1)*interv;
            end
            %% accuracy estimation for single-label prediction 
            Yones=zeros(length(X{l}),nm); Yones(:,l)=1;
            for i = 1:model.nThreshs,
                idxDet = Wc{l} <= thresholds(l,i);
                TP=sum(sum(idxDet & Yones));
                TN=sum(sum(~idxDet & ~Yones));
                FP=sum(sum(idxDet & ~Yones));
                FN=sum(sum(~idxDet & Yones));

                if model.accuracyglobal,
                    %%% global accuracy:
                    scoresA(l,i)= (TN+TP)./(TP+TN+FP+FN);
                else
                    %%% weighted accuracy (for imbalanced data sets). 
                    scoresA(l,i)= (TP/(TP+FN) + TN/(TN+FP))/2;
                end

                scoresP(l,i) = TP./(TP+FP);
                scoresR(l,i) = TP./(TP+FN);

                if isnan(scoresP(l,i)), scoresP(l,i) = 0; end
                if isnan(scoresR(l,i)), scoresR(l,i) = 0; end
                if isnan(scoresA(l,i)), scoresA(l,i) = 0; end
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
        display(sprintf('Classifying %d sequences in %d gesture classes ...',length(Y.L),nm));
        Wc = zeros(length(Y.L),nm);
    %     KT1 = []; KT2 = [];
        for s = 1:length(Y.L)      % save each sequence vs model dtw costs
            if s < length(Y.L)
                seq = X(Y.seg(s):Y.seg(s+1)-1,:);
            else
                seq = X(Y.seg(s):Y.seg(s+1),:);
            end
            if ~isempty(model.D)
                % 1:high-arm-wave , 17: side-kick
    %             if isempty(KT1) && strcmp(Yc.cnames(s),'high-arm-wave'), KT1 = getUpdatedCosts(seq,model.SM); end;
    %             if isempty(KT2) && strcmp(Yc.cnames(s),'side-kick'), KT2 = getUpdatedCosts(seq,model.SM); end;
    %             if ~isempty(KT1) && ~isempty(KT2)
    %                 display();
    %             end
                if model.darwin
                    KT = model.KT;
                else
                    KT = getUpdatedCosts(seq,model.SM);
                end
            end
            for k = 1:nm
                %% Compute the costs of the test sequence in terms of SM
                M = model.M{k};
                if ~isempty(model.D)
                    if ~iscell(M)
                        W = single(dtwc(seq,M,true,Inf,model.D,model.KM{k},KT));
                    else
                        if ~isempty(M{model.k})
                            W = single(dtwc(seq,M{model.k},true,Inf,model.D,model.KM{k},KT));
                        end
                    end
                else
                    if ~iscell(M)
                        if ~model.pdtw,
                            W = single(dtwc(seq,M,true));
                        else
                            Pql=zeros(size(seq,1),size(M,1));
                            for hh=1:size(M,1),
                                if ~isempty(model.lmodel(k,hh).obj),
                                    Pql(:,hh)= mixGaussLogprob(model.lmodel(k,hh).obj,seq);
        %                             Pql(:,hh)= log(pdf(model.lmodel(k,hh).obj,seq));
                                    Pql(isinf(Pql(:,hh)))=0;
                                end
                            end
                            noze = sum(Pql)~=0;
                            Pql(Pql==0) = mean(mean(Pql(:,noze)));
                            maval = max(max(abs(Pql)));
        %                     Dima  = (1-Pql)./maval;
        %                     DD=pdist2(seq,M);
        %                     DD = DD./max(max(DD));
                            Dima = (1-Pql)./maval.*pdist2(seq,M);
                            W = single(dtw3(seq,M,true,Inf,Dima));
                        end
                    else
                        if ~isempty(M{model.k})
                            W = single(dtwc(seq,M{model.k},true));
                        end
                    end
                end
                Wc(s,k) = W(end,end);
            end
        end
        if isempty(model.bestThs)
            for k = 1:nm     % save dtw costs of each non-iddle sequence vs its model
                interv = (max(Wc(Y.L==k,k))-min(Wc(Y.L==k,k)))/model.nThreshs;
                tMin = min(Wc(Y.L==k,k));
                if interv == 0, interv = tMin*2/model.nThreshs; end
                thresholds(k,:) = tMin + ((1:model.nThreshs)-1)*interv;
            end
        end
    %     scoresP = zeros(length(Y.L),model.nThreshs); scoresR = zeros(length(Y.L),model.nThreshs); scoresA = zeros(length(Y.L),model.nThreshs);
        %% accuracy estimation for single-label prediction 
        Yones=zeros(size(Y.L,2),nm);
        for i=1:nm,
            Yones(Y.L==i,i)=1;
        end
        for i = 1:model.nThreshs,
            for s = 1:nm,            
                idxDet = Wc(:,s)<=thresholds(s,i);  
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
    %     for s = 1:length(Y.L)
    %         for i = 1:model.nThreshs
    %             TP = 0; FP = 0; FN = 0; TN = 0;  % NO SE CONSIDERAN LOS TN, DECIDIR SI LOS USAMOS SEGUN LO HABLADO
    %             idxDet = Wc(s,:) < thresholds(:,i)';
    %             if Y.L(s) < nm+1    % iddle gesture is ignored if it wasn't learnt
    %                 if idxDet(Y.L(s))
    %                     TP = TP + 1;    % DECIDISION DE LA ASIGNACION DE VERDADEROS POSITIVOS
    %                     if sum(idxDet) > 1, FP = FP + sum(idxDet)-1; end
    %                 else
    %                     FN = FN + 1; FP = FP + sum(idxDet); 
    %                 end
    %             else
    %                 FP = FP + sum(idxDet); 
    %             end
    % %             if sum(~idxDet)       % DECISION PARA ASIGNACION DE VERDADEROS NEGATIVOS
    % %                 if ~idxDet(Y.L(s))
    % %                     TN = TN + sum(~idxDet)-1;
    % %                 else
    % %                     TN = TN + sum(~idxDet);
    % %                 end
    % %             end
    %             scoresP(s,i) = TP./(TP+FP);
    %             scoresR(s,i) = TP./(TP+FN);
    %             scoresA(s,i) = (TP)./(TP+FN+FP);    %%% (TP/(TP+FN) + TN/(TN+FP))/2   % METRICA PARA EL BALANCEO ENTRE CLASES
    %             if isnan(scoresP(s,i)), scoresP(s,i) = 0; end
    %             if isnan(scoresR(s,i)), scoresR(s,i) = 0; end
    %             if isnan(scoresA(s,i)), scoresA(s,i) = 0; end
    %         end
    %     end
        [bestScores(1),bestThsPos(1)] = max(mean(scoresP));
        [bestScores(2),bestThsPos(2)] = max(mean(scoresR));
        [bestScores(3),bestThsPos(3)] = max(mean(scoresA));
        score = bestScores(model.score2optim);
        model.bestThs = zeros(1,k); model.nThreshs = 1;
        for i = 1:k
            model.bestThs(i) = thresholds(i,bestThsPos(model.score2optim));
        end
    end
else
    %% Spotting
    display(sprintf('Spotting the %d classes in a sequence of %d frames ...',nm,length(Y.Lfr)));
    global DATATYPE; global NAT;
    %% Compute the costs of the test sequence in terms of SM 
    if ~isempty(model.D)
        display('Computing the costs of the test sequence in terms of SM ...');
        if ~model.darwin
            model.KT = getUpdatedCosts(X,model.SM);
        end
    end
    scoresO = zeros(nm,model.nThreshs);
    prexp = false(model.nThreshs,length(Y.Lfr),nm);
    bestScores = zeros(nm,4); bestThsPos = zeros(nm,4);
    %% Estimate scores & predictions of aligning X to each model
    for k = 1:nm
        GTtestk = Y.L == k; GTtestkFr = Y.Lfr == k;
        if ~any(GTtestkFr)
            warning('g:missedLabel','Label %d is missing in the test sequence',k);
        end
        W = [];
        if ~isempty(model.D)
            if ~iscell(model.M{k})
                W = single(dtwc(X,model.M{k},false,Inf,model.D,model.KM{k},model.KT));
            else
                if ~isempty(model.M{k}{model.k})
                    W = single(dtwc(X,model.M{k}{model.k},false,Inf,model.D,model.KM{k},model.KT));
                end
            end
        else
            if ~iscell(model.M{k})
                if ~model.pdtw,
                    W = single(dtwc(X,model.M{k},false));
                else
                    Pql=zeros(size(X,1),size(model.M{k},1));
                    for hh=1:size(model.M{k},1),
                        if ~isempty(model.lmodel(k,hh).obj),
                            Pql(:,hh)= mixGaussLogprob(model.lmodel(k,hh).obj,X);
%                             Pql(:,hh)= log(pdf(model.lmodel(k,hh).obj,X));
                            Pql(isinf(Pql(:,hh)))=0;
                        end
                    end
                    noze= sum(Pql)~=0;
                    Pql(Pql==0)=mean(mean(Pql(:,noze)));
                    maval=max(max(abs(Pql)));
%                     Dima  = (1-Pql)./maval;
%                     DD=pdist2(X,model.M{k});
%                     DD = DD./max(max(DD));
                    Dima  = (1-Pql)./maval.*pdist2(X,model.M{k});
                    W=single(dtw3(X,model.M{k},false,Inf,Dima));
                end
            else
                if ~isempty(model.M{k}{model.k})
                    W = single(dtwc(X,model.M{k}{model.k},false));
                end
            end
        end
        if isempty(model.bestThs)
            interv = (max(W(end,2:end))-min(W(end,2:end)))/model.nThreshs; 
            tMin = min(W(end,2:end));
            if interv == 0, interv = tMin*2/model.nThreshs; end
            thresholds(k,:) = tMin + ((1:model.nThreshs)-1)*interv;
        end
        if ~isempty(W)
%             detSeqLog3 = false(1,length(X));
            idxEval = [];            
            for i = 1:model.nThreshs                
                idx = find(W(end,:) <= thresholds(k,i));
                idx(ismember(idx,idxEval)) = [];
                %% Old, much slower
%                 tic;
%                 toc;
%                 tic;
%                 for j = 1:length(idx)
%     %                 if detSeqLog(idx(j)-1)==0
%         %                 fprintf('%d ',idx(j)-1);
%     %                     fprintf('testing with threshold %.2f\n',idx(j));
%     %                     [in,fi,~] = aligngesture([],W(:,1:idx(j)));
%                         [in,fi] = detectSeqC(W(:,1:idx(j)));
%     %                     [~,in3,fi3] = getDTWcseq(W(:,1:idx(j)));
%                         if ~isequal(in,in2,in3) || ~isequal(fi,fi2,fi3)
%                             disp('');
%                         end
%                         if length(in) > 1 || length(fi) > 1
%                             error('Start and end of the gesture must be scalars');
%                         end
%                         detSeqLog3(in:fi) = 1;
%     %                 end
%                 end
%                 toc;
                %% This is much faster
                prexp(i,:,k) = getDetectedSeqs_c(W,int32(idx),prexp(i,:,k),model.maxWlen);
                %%
                % to compensate for the offset of deep-features
                if strcmp(DATATYPE,'chalearn2014') && NAT == 3
                    prexp(i,:,k)=([prexp(i,6:end,k),0,0,0,0,0]);    % correct offset of deep features
                end
                idxEval = unique([idxEval idx(prexp(i,idx-1,k)==true)]);
%                 if ~isequal(detSeqLog3,detSeqLog)
%                     find(detSeqLog3~=detSeqLog)
%                     if sum(detSeqLog3~=detSeqLog) > 1
%                         error();
%                     end
%                 end
                detSw = prexp(i,:,k);
                scoresO(k,i) = sum(GTtestkFr & detSw)./sum(GTtestkFr | detSw);     % overlap (Jaccard Index)
                detSw = getActivations(detSw, GTtestkFr, Y.seg, model);
                
                % only for MADX database (recognition)
                if strcmp(DATATYPE,'mad1') || strcmp(DATATYPE,'mad2') ...
                        || strcmp(DATATYPE,'mad3') || strcmp(DATATYPE,'mad4') ...
                        || strcmp(DATATYPE,'mad5') 
                    [~,~,R] = estimate_overlap_madold(GTtestk, detSw, model.minOverlap);
                    scoresP(k,i) = R.prec2;  % Precision
                    scoresR(k,i) = R.rec2;    % Recall
                else
                    scoresP(k,i) = sum(GTtestk & detSw)./sum(GTtestk & detSw | ~GTtestk & detSw);  % Precision
                    scoresR(k,i) = sum(GTtestk & detSw)./sum(GTtestk & detSw | GTtestk & ~detSw);  % Recall
                end
                scoresA(k,i) = sum(GTtestk & detSw | ~GTtestk & ~detSw)./sum(GTtestk & detSw | GTtestk & ~detSw | ~GTtestk & detSw | ~GTtestk & ~detSw);   % Accuracy
                if isnan(scoresO(k,i)), scoresO(k,i) = 0; end                
                if isnan(scoresP(k,i)), scoresP(k,i) = 0; end
                if isnan(scoresR(k,i)), scoresR(k,i) = 0; end
                if isnan(scoresA(k,i)), scoresA(k,i) = 0; end
            end
        end
        [bestScores(k,1),bestThsPos(k,1)] = max(scoresO(k,:));
        [bestScores(k,2),bestThsPos(k,2)] = max(scoresP(k,:));
        [bestScores(k,3),bestThsPos(k,3)] = max(scoresR(k,:));
        [bestScores(k,4),bestThsPos(k,4)] = max(scoresA(k,:));        
        %% save mean scores and learnt thresholds (allow multi-label assignment per frame)
%         score = mean(bestScores(:,model.score2optim));
%         bestScores = mean(bestScores)
        
%         gtF = Y.Lfr;
%         save('predictions.mat','predictionsF','gtF');
%         imagesc(confusionmat(predictionsF,gtF)); colormap(hot);
    
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
%         mean(bestScores)
        clear bestScores;
        
        %%% if we want to plot predictions vs GT
        % plot(predictions);
        % gcf;hold on;
        % plot(nY,':r');

        % now estimate overlap, prec, rec, and f1 (no accuracy right now)
        [~,~,R,ovlp] = estimate_overlap_mad(Y,predictions,model.minOverlap);
        % Assign Overlap, Precision, Recall, F1-score
        bestScores(1) = ovlp;  bestScores(2) = R.prec; bestScores(3) = R.rec; bestScores(4) = (2.*R.rec.*R.prec)./(R.rec + R.prec);
%         bestScores
        for j=1:length(bestScores), if isnan(bestScores(j)), bestScores(j)=0; end; end
        
        switch model.score2optim
            case 1, score=ovlp;
            case 2, score=R.prec;
            case 3, score=R.rec;
            case 4, score=(2.*R.rec.*R.prec)./(R.rec + R.prec); % F-1 score
        end
        clear prexp R ovlp ofinscores ofin cupred cues;
        if isempty(model.bestThs)
            model.bestThs = zeros(1,k); model.nThreshs = 1; 
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