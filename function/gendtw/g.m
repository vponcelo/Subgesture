function [model,score,bestScores,predictions] = g(model,Xc,Yc)
% Output:
%   model
%   score(default): mean overlap for the k detected gestures
% Input:
%   model: model parameters
%   Xc: data to test against (cell or whole sequence)
%   Yc: labels of X (cell or whole sequence)

global DATATYPE; global NAT;

if iscell(Xc)
    error('g:classErr','Input data cannot be a cell here, need to specify this functionality');
end
sw = model.sw;
if sw > length(Yc.Lfr)
    error('g:swErr','sliding window is longer than the validation sequence');
end
ns = 1;     % number of sliding windows, more than one ns > 1 is not recommended
if sw == 0
    sw = length(Xc)-1;
    r = 1;
else
    sw = model.sw;
    r = inf;
    while any(r > length(Yc.Lfr)-sw)
        r = randperm(round(length(Yc.Lfr)),ns);
    end
end
  
thresholds = cell(1,ns);
scoresO = cell(1,ns);
if model.classification
    scoresP = cell(1,ns); scoresR = cell(1,ns); scoresA = cell(1,ns);
else
    predictions = cell(1,ns);
end
nm = length(model.M);

for isw = 1:ns  % sliding window (ns > 1 not recommended)
    % last sliding window is the size of the sequence
    predictions{isw} = [];
    seg=r(isw):min(r(isw)+sw,length(Yc.Lfr));
    X=Xc(seg,:);
    Y.Lfr=Yc.Lfr(seg);
    if model.classification
        Y.L=Y.Lfr; d=diff(Y.Lfr); Y.L(d==0)=[]; Y.seg=[1 find(d~=0)+1 length(Y.Lfr)];
    end
    display(sprintf('Evaluating %d of %d sequence of %d frames ...',isw,ns,length(Y.Lfr)));
    %% Compute the costs of the test sequence in terms of SM 
    if ~isempty(model.D)
        display('Computing the costs of the test sequence in terms of SM ...');
%         tic;
        model.KT = getUpdatedCosts(X,model.SM);
%         toc;        
    end
        
    %% Learn threshold cost parameters for each gesture
    thresholds{isw} = cell(1,nm);
    scoresO{isw} = cell(1,nm);
    if model.classification
        scoresP{isw} = cell(1,nm); scoresR{isw} = cell(1,nm); scoresA{isw} = cell(1,nm);
    end
    
    for k = 1:nm
        if model.classification, GTtestk = Y.L == k; end
        GTtestkFr = Y.Lfr == k;
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
        swthreshs = zeros(1,model.nThreshs); swOvs = zeros(1,model.nThreshs);
        if model.classification
            swPrecs = zeros(1,model.nThreshs); swRecs = zeros(1,model.nThreshs);
            swAccs = zeros(1,model.nThreshs);
        end
        detSeqLog = false(model.nThreshs,length(GTtestkFr));
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
                if strcmp(DATATYPE,'chalearn2014') && NAT == 3
                    detSeqLog(i,:)=([detSeqLog(i,6:end),0,0,0,0,0]);    % correct offset of deep features
                end
                idxEval = unique([idxEval idx(detSeqLog(i,idx-1)==true)]);
%                 if ~isequal(detSeqLog3,detSeqLog)
%                     find(detSeqLog3~=detSeqLog)
%                     if sum(detSeqLog3~=detSeqLog) > 1
%                         error();
%                     end
%                 end
                swOvs(i) = sum(GTtestkFr & detSeqLog(i,:))./sum(GTtestkFr | detSeqLog(i,:));     % overlap (Jaccard Index)
                if isnan(swOvs(i)), swOvs(i) = 0; end
                if model.classification
                    detSw = getActivations(detSeqLog(i,:), GTtestkFr, Y.seg, model.minOverlap);
                    swPrecs(i) = sum(GTtestk & detSw)./sum(GTtestk & detSw | ~GTtestk & detSw);  % Precision
                    if isnan(swPrecs(i)), swPrecs(i) = 0; end
                    swRecs(i) = sum(GTtestk & detSw)./sum(GTtestk & detSw | GTtestk & ~detSw);  % Recall
                    if isnan(swRecs(i)), swRecs(i) = 0; end
                    swAccs(i) = sum(GTtestk & detSw | ~GTtestk & ~detSw)./sum(GTtestk & detSw | GTtestk & ~detSw | ~GTtestk & detSw | ~GTtestk & ~detSw);   % Accuracy
                    if isnan(swAccs(i)), swAccs(i) = 0; end
                end
            end
        end
        thresholds{isw}{k} = swthreshs;
        scoresO{isw}{k} = swOvs;
        if model.classification
            scoresP{isw}{k} = swPrecs; scoresR{isw}{k} = swRecs; scoresA{isw}{k} = swAccs;
        else
            [~,pos] = max(scoresO{isw}{k});
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

if ~model.classification,    
    bestScores = zeros(k,ns,1); bestThsPos = zeros(k,ns,1); 
else
    bestScores = zeros(k,ns,4); bestThsPos = zeros(k,ns,4);
end
for i = 1:k
    for j = 1:ns
        [bestScores(i,j,1),bestThsPos(i,j)] = max(scoresO{j}{i});
        if model.classification
            [bestScores(i,j,2),bestThsPos(i,j,2)] = max(scoresP{j}{i});
            [bestScores(i,j,3),bestThsPos(i,j,3)] = max(scoresR{j}{i});
            [bestScores(i,j,4),bestThsPos(i,j,4)] = max(scoresA{j}{i});            
        end    
    %     gtF = Y.Lfr;
    %     save('predictions.mat','predictionsF','gtF');
    %     imagesc(confusionmat(predictionsF,gtF)); colormap(hot);
    end
end
if ~model.classification
    score = mean(bestScores(:,:,1));
else
    score = mean(bestScores(:,:,model.score2Optim));     
end
score = max(score);     % not reliable if ns > 1
if ns == 1
    bestScores = reshape(mean(bestScores(:,ns,:)),[1 size(bestScores,3)]);
end
% try
%     seg=r(p):min(r(p)+sw,length(Yc.Lfr));
%     Y.Lfr=Yc.Lfr(seg);
%     plotmistakes(predictions{p},Y,1);
% catch e
%     display(e.message);
% end
% close gcf;
    
if strcmp(model.scoreMeasure,'overlap')
    model.bestThs = zeros(1,k);
    for i = 1:k
        ths = zeros(1,ns);
        for j = 1:ns
            ths(j) = thresholds{j}{i}(bestThsPos(i,j,model.score2Optim));
        end
        model.bestThs(i) = mean(ths);   % not reliable if ns > 1
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