function [model,score,predictions] = g(model,Xc,Yc)
% Output:
%   model
%   score(default): mean overlap for the k detected gestures
% Input:
%   model: model parameters
%   Xc: data to test against (cell or whole sequence)
%   Yc: labels of X (cell or whole sequence)

if iscell(Xc)
    error('g:classErr','Input data cannot be a cell here, need to specify this functionality');
end
sw = model.sw;
if sw > length(Yc.Lfr)
    error('g:swErr','sliding window is longer than the validation sequence');
end
if sw == 0
    sw = length(Xc)-1;
    ns = 1;
    r=1;
else
    sw = model.sw;
    ns = 1;     %round(length(Yc.Lfr)/sw)     % changing this value to 1 evaluates the sequence once
    r = inf;
    while any(r > length(Yc.Lfr)-sw)
        r=randperm(round(length(Yc.Lfr)),ns);
    end
end
  
thresholds = cell(1,ns);
scores = cell(1,ns);
if ~model.classification && strcmp(model.scoreMeasure,'levenshtein')
    predictions = cell(1,ns);
end

for isw = 1:ns  % sliding window
    % last sliding window is the size of the sequence
    display(sprintf('Evaluating %d of %d sequences of %d frames ...',isw,ns,sw));
    predictions{isw} = [];
    seg=r(isw):min(r(isw)+sw,length(Yc.Lfr));
    X=Xc(seg,:);
    Y.Lfr=Yc.Lfr(seg);
    if model.classification
        Y.L=Yc.Lfr; d=diff(Yc.Lfr); Y.L(d==0)=[]; Y.seg=[1 find(d~=0)+1];
    end
    
    %% Compute the costs of the test sequence in terms of SM 
    if ~isempty(model.D)
        display('Computing the costs of the test sequence in terms of SM ...');
%         tic;
        model.KT = getUpdatedCosts(X,model.SM);
%         toc;        
    end
        
    %% Learn threshold cost parameters for each gesture
    thresholds{isw} = cell(1,length(model.M));
    scores{isw} = cell(1,length(model.M));
    
    for k = 1:length(model.M)
        thresholds{isw}{k} = [];
        scores{isw}{k} = [];
    end
    
    for k = 1:length(model.M)
        
    %     if ~isempty(model.D)    
    %         display('Computing the partial costs of the test sequence in terms of SM ...');
    %         model.KT = getUpdatedCosts(X,model.SM);
    %     end
        %         display(sprintf('Testing threshold cost parameters for gesture %d ...',k));
        if model.classification
            GTtestk = Y.L == k;
        else
            GTtestk = Y.Lfr == k;
        end
        if ~any(GTtestk)
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
            TOL_THRESH = 0.001;
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
            TOL_THRESH = 0.01;
        end
        if ~isempty(model.bestThs)
            model.nThreshs = 1;        
        else
            if ~isempty(W)
                interv = (max(W(end,2:end))-min(W(end,2:end)))/model.nThreshs;
                while interv < TOL_THRESH && model.nThreshs > 1
                    model.nThreshs = round(model.nThreshs/2);
                    interv = (max(W(end,2:end))-min(W(end,2:end)))/model.nThreshs;
                    display(sprintf('Decreasing number of test thresholds to %d',model.nThreshs));
                end
            end
        end
        swthreshs = zeros(1,model.nThreshs);
        swScores = zeros(1,model.nThreshs);
        detSeqLog = false(model.nThreshs,length(Y.Lfr));
        if ~isempty(W)
            tMin = min(W(end,2:end));
%             detSeqLog3 = false(1,length(X));
            idxEval = [];
            for i = 1:model.nThreshs
                if isempty(model.bestThs)
                    swthreshs(i) = tMin + (i-1)*interv;
                else
                    swthreshs(i) = model.bestThs(k);
                end
                idx = find(W(end,:) <= swthreshs(i));
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
                detSeqLog(i,:) = getDetectedSeqs_c(W,int32(idx),detSeqLog(i,:),model.maxWlen);
                %%
                % to compensate for the offset of deep-features
                detSeqLog(i,:)=([detSeqLog(i,6:end),0,0,0,0,0]);
                idxEval = unique([idxEval idx(detSeqLog(i,idx-1)==true)]);
%                 if ~isequal(detSeqLog3,detSeqLog)
%                     find(detSeqLog3~=detSeqLog)
%                     if sum(detSeqLog3~=detSeqLog) > 1
%                         error();
%                     end
%                 end

                if model.classification
                    d = diff(detSeqLog(i,:)); idxL = find(d~=0); 
                    detSw = zeros(1,length(GTtestk));
                    for l = 1:length(idxL)
                        cLabel = sum(Y.seg < idxL(l));
                        if cLabel > 0
                            detSeqLog(i,cLabel) = 1;
                        end
                    end
%                     swScores(i) = sum(GTtestk & detSw)./sum(GTtestk & detSw | ~GTtestk & detSw);  % Precision
%                     swScores(i) = sum(GTtestk & detSw)./sum(GTtestk & detSw | GTtestk & ~detSw);  % Recall
                    swScores(i) = sum(GTtestk & detSw | ~GTtestk & ~detSw)./sum(GTtestk & detSw | GTtestk & ~detSw | ~GTtestk & detSw | ~GTtestk & ~detSw);   % Accuracy
                else
                    swScores(i) = sum(GTtestk & detSeqLog(i,:))./sum(GTtestk | detSeqLog(i,:));     % overlap (Jaccard Index)
                end
            end
        end
        thresholds{isw}{k} = swthreshs;
        scores{isw}{k} = swScores;
        if ~model.classification && strcmp(model.scoreMeasure,'levenshtein')
            [~,pos] = max(scores{isw}{k});
            idx = find(detSeqLog(pos,:) == 1);
            if ~isempty(idx)
                inF = idx(1);
                for i = 2:length(idx)
                    if idx(i)-idx(i-1) > 1
                        endF = idx(i-1);
                        predictions{isw} = [predictions{isw} inF endF k];
                        inF = idx(i);
                    end
                    if i == length(idx)
                        endF = idx(i);
                        predictions{isw} = [predictions{isw} inF endF k];
                        % Last ones when inF == endF will be also added.
                    end
                end
            else
                predictions{isw} = 0;
            end
        end
    end
end

bestScores = cell(1,ns);
score = zeros(1,ns);
for j = 1:ns
    bestScores{j} = zeros(1,k);
    for i = 1:k        
        [bestScores{j}(i),~] = max(scores{j}{i});
    end
    score(j) = mean(bestScores{j});
    %     gtF = Y.Lfr;
    %     save('predictions.mat','predictionsF','gtF');
    %     imagesc(confusionmat(predictionsF,gtF)); colormap(hot);
end
[score,p] = max(score); % should be the mean

% try
%     seg=r(p):min(r(p)+sw,length(Yc.Lfr));
%     Y.Lfr=Yc.Lfr(seg);
%     plotmistakes(predictions{p},Y,1);
% catch e
%     display(e.message);
% end
% close gcf;
    
if strcmp(model.scoreMeasure,'overlap')
    model.bestThs = zeros(1,length(thresholds{p}));
    for i = 1:k
        [~,pos] = max(scores{p}{i});
        model.bestThs(i) = thresholds{p}{i}(pos);
    end
%     gtF = Y.Lfr;
%     save('predictions.mat','predictionsF','gtF');
%     imagesc(confusionmat(predictionsF,gtF)); colormap(hot);
elseif strcmp(model.scoreMeasure,'levenshtein')
    predLabels = []; lints = [];
    for i = 1:size(X,1)
        l = find(predictions{p}(1:3:end) == i);
        while ~isempty(l) && ~any(ismember(l,lints))
            [minFi,pos] = min(predictions{p}(3*(l-1)+2));
            predLabels = [predLabels predictions{p}(3*l(pos))];
            l(pos) = [];
            lint = find(predictions{p}(1:3:end) > i & predictions{p}(2:3:end) < minFi);
            lint(ismember(lint,lints)) = [];
            lints = [lints lint];
            while ~isempty(lint)
                [~,posInt] = min(predictions{p}(3*(lint-1)+1));
                predLabels = [predLabels predictions{p}(3*lint(posInt))];
                lint(posInt) = [];
            end
        end
    end
    score = levenshtein(predLabels,Y.L);        
else
    error('g:errorMeasure','Invalid measure. Check scoreMeasure parameter');
end