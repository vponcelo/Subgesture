function [model,score,predictions] = g(model,Xc,Yc)
% Output:
%   model
%   score(default): mean overlap for the k detected gestures
% Input:
%   model: model parameters
%   X: data to test against
%   Y: labels of X

%% Compute the costs of the test sequence in terms of SM 
if ~iscell(Xc)
    X = Xc;
    Y = Yc;
end
if model.sw~=0,
    seglen=model.sw;
    r=randperm(round(length(Y.Lfr)./2));
    segment=r(1):min(r(1)+seglen,length(Y.Lfr));
    X=X(segment,:);
    Y.Lfr=Y.Lfr(segment);
end

if ~isempty(model.D)    
    display('Computing the costs of the test sequence in terms of SM ...');
    tic;
    model.KT = getUpdatedCosts(X,model.SM);
    toc;        
end
        
%% Learn threshold cost parameters for each gesture
folds = 1;
thresholds = cell(1,length(model.M));
overlaps = cell(1,length(model.M));
predictions = [];
predictionsF = zeros(1,length(X));
for k = 1:length(model.M)
    if iscell(Xc)
        X = Xc{k};
        Y = Yc{k};
    end
%     if ~isempty(model.D)    
%         display('Computing the partial costs of the test sequence in terms of SM ...');
%         model.KT = getUpdatedCosts(X,model.SM);
%     end
    fprintf('Testing threshold cost parameters for gesture %d ...\n',k);
    GTtestk = Y.Lfr == k;
    if ~any(GTtestk)
        error('g:missedLabel','Label %d is missing in the test sequence',k);
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
            end
        end
    end
    thresholds{k} = zeros(folds,model.nThreshs);
    overlaps{k} = zeros(folds,model.nThreshs);
    for K = 1:folds
        detSeqLog = false(model.nThreshs,length(X));
        if ~isempty(W)
            tMin = min(W(end,2:end));
%             detSeqLog3 = false(1,length(X));
            idxEval = [];
            for i = 1:model.nThreshs                
                if isempty(model.bestThs)
                    thresholds{k}(K,i) = tMin + (i-1)*interv;
                else
                    thresholds{k}(K,i) = model.bestThs(k);
                end
                idx = find(W(end,:) <= thresholds{k}(K,i));
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
                overlaps{k}(K,i) = sum(GTtestk & detSeqLog(i,:))/sum(GTtestk | detSeqLog(i,:));
            end
        end
        [~,pos] = max(overlaps{k});
        idx = find(detSeqLog(pos,:) == 1);
        predictionsF(idx) = k;
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
    end
end
display('Done!');
try
plotmistakes(predictions,Y,1);
catch
    lasterr
end
close gcf;

if strcmp(model.scoreMeasure,'overlap')
    bestOverlaps = zeros(1,length(overlaps));
    model.bestThs = zeros(1,length(thresholds));
    for i = 1:k
        [bestOverlaps(i),pos] = max(overlaps{i});
        model.bestThs(i) = thresholds{i}(pos);
    end
    model.bestThs
    score = mean(bestOverlaps);
%     gtF = Y.Lfr;
%     save('predictions.mat','predictionsF','gtF');
%     imagesc(confusionmat(predictionsF,gtF)); colormap(hot);
elseif strcmp(model.scoreMeasure,'levenshtein')
    predLabels = []; lints = [];
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
else
    error('g:errorMeasure','Invalid measure. Check scoreMeasure parameter');
end

% function [model,score,predictions] = g(model,X,Y)
% % Output:
% %   model
% %   score(default): mean overlap for the k detected gestures
% % Input:
% %   model: model parameters
% %   X: data to test against
% %   Y: labels of X
% 
% %sw = length(X);
% sw = model.sw;
% if sw == 0
%     sw = length(X);
%     ns = 1;
% else
%     sw = model.sw;
%     ns = 3;     % changing this value to 1 evaluates the sequence once
% end
% 
% thresholds = cell(1,ns);
% overlaps = cell(1,ns);
% predictions = cell(1,ns);
% 
% for isw = 1:ns  % sliding window
%     % last sliding window is the size of the sequence
%     display(sprintf('Evaluating %d of %d sequences of %d frames ...',isw,ns,sw));
%     predictions{isw} = [];
%     s = 1; e = length(X);
%     while e >= isw*length(X)/ns && model.sw > 0
%         s = randi([round((isw-1)*length(X)/ns+1) round(isw*length(X)/ns)]);
%         e = s + sw;
%     end
%     
%     %% Compute the costs of the test sequence in terms of SM 
%     if ~isempty(model.D)
%         display('Computing the costs of the test sequence in terms of SM ...');
%         tic;
%         model.KT = getUpdatedCosts(X(s:e,:),model.SM);
%         toc;
%     end
% 
%     %% Learn threshold cost parameters for each gesture
%     folds = 1;
%     thresholds{isw} = cell(1,length(model.M));
%     overlaps{isw} = cell(1,length(model.M));
%     
%     for k = 1:length(model.M)
%         thresholds{isw}{k} = [];
%         overlaps{isw}{k} = [];
%     end
%     
%     for k = 1:length(model.M)
%     %     if ~isempty(model.D)    
%     %         display('Computing the partial costs of the test sequence in terms of SM ...');
%     %         model.KT = getUpdatedCosts(X(s:e,:),model.SM);
%     %     end
%         display(sprintf('Testing threshold cost parameters for gesture %d ...',k));
%         GTtestk(s:e) = Y.Lfr(s:e) == k;
%     %     if ~any(GTtestk)
%     %         error('g:missedLabel','Label %d is missing in the test sequence',k);
%     %     end
%         W = [];
%         if ~isempty(model.D)  
%             if ~iscell(model.M{k})
%                 W = single(dtwc(X(s:e,:),model.M{k},false,Inf,model.D,model.KM{k},model.KT));
%             else
%                 if ~isempty(model.M{k}{model.k})
%                     W = single(dtwc(X(s:e,:),model.M{k}{model.k},false,Inf,model.D,model.KM{k},model.KT));
%                 end
%             end
%             TOL_THRESH = 0.001;
%         else
%             if ~iscell(model.M{k})
%                 W = single(dtwc(X(s:e,:),model.M{k},false));
%             else
%                 if ~isempty(model.M{k}{model.k})
%                     W = single(dtwc(X(s:e,:),model.M{k}{model.k},false));            
%                 end
%             end        
%             TOL_THRESH = 0.01;
%         end
%         if ~isempty(model.bestThs)
%             model.nThreshs = 1;
%         else
%             if ~isempty(W)
%                 interv = (max(W(end,2:end))-min(W(end,2:end)))/model.nThreshs;
% %                 while interv < TOL_THRESH && model.nThreshs > 1
% %                     model.nThreshs = round(model.nThreshs/2);
% %                     interv = (max(W(end,2:end))-min(W(end,2:end)))/model.nThreshs;
% %                     display(sprintf('Incresing number of thresholds to %d',model.nThreshs));
% %                 end
%             end
%         end
%         swthreshs = zeros(folds,model.nThreshs);
%         swoverlaps = zeros(folds,model.nThreshs);
% %         thresholds{k} = zeros(folds,model.nThreshs);
% %         overlaps{k} = zeros(folds,model.nThreshs);
%         for K = 1:folds
%             detSeqLog = false(model.nThreshs,length(GTtestk));
%             if ~isempty(W)
%                 tMin = min(W(end,2:end));
%     %             detSeqLog3 = false(1,length(X(s:sw,:)));
%                 idxEval = [];
%                 for i = 1:model.nThreshs
%                     if isempty(model.bestThs)
%                         swthreshs(K,i) = tMin + (i-1)*interv;
%                     else
%                         swthreshs(K,i) = model.bestThs(k);
%                     end
%                     idx = find(W(end,:) <= swthreshs(K,i));
%                     idx(ismember(idx,idxEval)) = [];
%                     detSeqLog(i,:) = getDetectedSeqs_c(W,int32(idx),detSeqLog(i,:),model.maxWlen);
%     %                 tic;
%     %                 toc;
%     %                 tic;
%     %                 for j = 1:length(idx)
%     %     %                 if detSeqLog(idx(j)-1)==0
%     %         %                 fprintf('%d ',idx(j)-1);
%     %     %                     fprintf('testing with threshold %.2f\n',idx(j));
%     %     %                     [in,fi,~] = aligngesture([],W(:,1:idx(j)));
%     %                         [in,fi] = detectSeqC(W(:,1:idx(j)));
%     %     %                     [~,in3,fi3] = getDTWcseq(W(:,1:idx(j)));
%     %                         if ~isequal(in,in2,in3) || ~isequal(fi,fi2,fi3)
%     %                             disp('');
%     %                         end
%     %                         if length(in) > 1 || length(fi) > 1
%     %                             error('Start and end of the gesture must be scalars');
%     %                         end
%     %                         detSeqLog3(in:fi) = 1;
%     %     %                 end
%     %                 end
%     %                 toc;
%                     idxEval = unique([idxEval idx(detSeqLog(i,idx-1)==true)]);            
%     %                 if ~isequal(detSeqLog3,detSeqLog)
%     %                     find(detSeqLog3~=detSeqLog)
%     %                     if sum(detSeqLog3~=detSeqLog) > 1
%     %                         error();
%     %                     end
%     %                 end
%                     %%%%% to compensate for the offset of deep-features
%                     detSeqLog(i,:)=([detSeqLog(i,6:end),0,0,0,0,0]);
%                     swoverlaps(K,i) = sum(GTtestk & detSeqLog(i,:))/sum(GTtestk | detSeqLog(i,:));
% %                     overlaps{k}(K,i) = sum(GTtestk & detSeqLog(i,:))/sum(GTtestk | detSeqLog(i,:));
%                 end
%             end
%             thresholds{isw}{k} = swthreshs(K,:);
%             overlaps{isw}{k} = swoverlaps(K,:);
%             [~,pos] = max(overlaps{isw}{k});
%             idx = find(detSeqLog(pos,:) == 1);
%             if ~isempty(idx)
%                 inF = idx(1);
%                 for i = 2:length(idx)
%                     if idx(i)-idx(i-1) > 1
%                         endF = idx(i-1);
%                         predictions{isw} = [predictions{isw} inF endF k];
%                         inF = idx(i);
%                     end
%                     if i == length(idx)
%                         endF = idx(i);
%                         predictions{isw} = [predictions{isw} inF endF k];
%                         % Last ones when inF == endF will be also added.
%                     end
%                 end
%             else
%                 predictions{isw} = 0;
%             end
%         end
%     end
% end
% % display('Done!');
% 
% bestOverlaps = cell(1,ns);
% score = zeros(1,ns);
% for j = 1:ns
%     bestOverlaps{j} = zeros(1,k);
%     for i = 1:k        
%         [bestOverlaps{j}(i),~] = max(overlaps{j}{i});
%     end
%     score(j) = mean(bestOverlaps{j});
%     %     gtF = Y.Lfr;
%     %     save('predictions.mat','predictionsF','gtF');
%     %     imagesc(confusionmat(predictionsF,gtF)); colormap(hot);
% end
% [score,p] = max(score);
% 
% if strcmp(model.scoreMeasure,'overlap')    
%     model.bestThs = zeros(1,length(thresholds{p}));
%     for i = 1:k
%         [~,pos] = max(overlaps{p}{i});
%         model.bestThs(i) = thresholds{p}{i}(pos);
%     end    
% elseif strcmp(model.scoreMeasure,'levenshtein')
%     predLabels = []; lints = [];
%     for i = 1:size(X,1)
%         l = find(predictions{p}(1:3:end) == i);
%         while ~isempty(l) && ~any(ismember(l,lints))
%             [minFi,pos] = min(predictions{p}(3*(l-1)+2));
%             predLabels = [predLabels predictions{p}(3*l(pos))];
%             l(pos) = [];
%             lint = find(predictions{p}(1:3:end) > i & predictions{p}(2:3:end) < minFi);
%             lint(ismember(lint,lints)) = [];
%             lints = [lints lint];
%             while ~isempty(lint)
%                 [~,posInt] = min(predictions{p}(3*(lint-1)+1));
%                 predLabels = [predLabels predictions{p}(3*lint(posInt))];
%                 lint(posInt) = [];
%             end
%         end
%     end
%     score = levenshtein(predLabels,Y.L);
% else
%     error('g:errorMeasure','Invalid measure. Check scoreMeasure parameter');
% end
