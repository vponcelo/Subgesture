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
if sw == 0
    sw = length(Xc)-1;
    ns = 1;
else
    sw = model.sw;
    ns = 3;     % changing this value to 1 evaluates the sequence once    
end

r = inf;
while any(r > length(Yc.Lfr)-sw)
    r=randperm(round(length(Yc.Lfr)),ns);
end

% X=X(seg(isw),:);
% Y.Lfr=Y.Lfr(seg(isw));
    
thresholds = cell(1,ns);
overlaps = cell(1,ns);
predictions = cell(1,ns);

for isw = 1:ns  % sliding window
    % last sliding window is the size of the sequence
    display(sprintf('Evaluating %d of %d sequences of %d frames ...',isw,ns,sw));
    predictions{isw} = [];
    seg=r(isw):min(r(isw)+sw,length(Yc.Lfr));
    X=Xc(seg,:);
    Y.Lfr=Yc.Lfr(seg);
    
    %% Compute the costs of the test sequence in terms of SM 
    if ~isempty(model.D)    
        display('Computing the costs of the test sequence in terms of SM ...');
%         tic;
        model.KT = getUpdatedCosts(X,model.SM);
%         toc;        
    end
        
    %% Learn threshold cost parameters for each gesture
    folds = 1;
    thresholds{isw} = cell(1,length(model.M));
    overlaps{isw} = cell(1,length(model.M));
    
    for k = 1:length(model.M)
        thresholds{isw}{k} = [];
        overlaps{isw}{k} = [];
    end
    
    for k = 1:length(model.M)
        
    %     if ~isempty(model.D)    
    %         display('Computing the partial costs of the test sequence in terms of SM ...');
    %         model.KT = getUpdatedCosts(X,model.SM);
    %     end
        %         display(sprintf('Testing threshold cost parameters for gesture %d ...',k));
        GTtestk = Y.Lfr == k;
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
                W = single(dtwc(X,model.M{k},false));
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
                    display(sprintf('Incresing number of test thresholds to %d',model.nThreshs));
                end
            end
        end
        swthreshs = zeros(folds,model.nThreshs);
        swoverlaps = zeros(folds,model.nThreshs);
        for K = 1:folds
            detSeqLog = false(model.nThreshs,length(GTtestk));
            if ~isempty(W)
                tMin = min(W(end,2:end));
    %             detSeqLog3 = false(1,length(X));
                idxEval = [];
                for i = 1:model.nThreshs                
                    if isempty(model.bestThs)
                        swthreshs(K,i) = tMin + (i-1)*interv;
                    else
                        swthreshs(K,i) = model.bestThs(k);
                    end
                    idx = find(W(end,:) <= swthreshs(K,i));
                    idx(ismember(idx,idxEval)) = [];
                    detSeqLog(i,:) = getDetectedSeqs_c(W,int32(idx),detSeqLog(i,:),model.maxWlen);
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
                    idxEval = unique([idxEval idx(detSeqLog(i,idx-1)==true)]);            
    %                 if ~isequal(detSeqLog3,detSeqLog)
    %                     find(detSeqLog3~=detSeqLog)
    %                     if sum(detSeqLog3~=detSeqLog) > 1
    %                         error();
    %                     end
    %                 end
                    %%%%% to compensate for the offset of deep-features
                    detSeqLog(i,:)=([detSeqLog(i,6:end),0,0,0,0,0]);
                    swoverlaps(K,i) = sum(GTtestk & detSeqLog(i,:))/sum(GTtestk | detSeqLog(i,:));
                end
            end
            thresholds{isw}{k} = swthreshs(K,:);
            overlaps{isw}{k} = swoverlaps(K,:);
            [~,pos] = max(overlaps{isw}{k});
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

bestOverlaps = cell(1,ns);
score = zeros(1,ns);
for j = 1:ns
    bestOverlaps{j} = zeros(1,k);
    for i = 1:k        
        [bestOverlaps{j}(i),~] = max(overlaps{j}{i});
    end
    score(j) = mean(bestOverlaps{j});
    %     gtF = Y.Lfr;
    %     save('predictions.mat','predictionsF','gtF');
    %     imagesc(confusionmat(predictionsF,gtF)); colormap(hot);
end
[score,p] = max(score);

% try
%     seg=r(p):min(r(p)+sw,length(Yc.Lfr));
%     Y.Lfr=Yc.Lfr(seg);
%     plotmistakes(predictions{p},Y,1);
% catch e
%     display(e.message);
% end
close gcf;
    
if strcmp(model.scoreMeasure,'overlap')
    model.bestThs = zeros(1,length(thresholds{p}));
    for i = 1:k
        [~,pos] = max(overlaps{p}{i});
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