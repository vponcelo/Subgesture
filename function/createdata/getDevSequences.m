function [seq,GT] = getDevSequences(X,Y,l,noise,secs,nSampGest)
% Obtain a frame sequence and labels given the whole batches
% Output: 
%   seq: test sequences generated for training and validation
%   GT: ground truths of the training and validation sequences
% Input:
%   X: data sequences
%   Y: data labels
%   l: list with the number of batches to consider
%   noise: flag indicating whether or not consider noise
%   secs: reference seconds for the test sequence
%   nSampGest: Number of samples per gesture for the test sequence   

seq = cell(1,2);
GT = cell(1,2);

for v = 1:length(seq)
    if isempty(l)
        seq{v} = zeros(Y{v}.seg(end),size(X{v},2));
        GT{v}.Lfr = zeros(1,Y{v}.seg(end));
        for j = 1:length(Y{v}.L)
            startSeq = Y{v}.seg(j);
            GT{v}.seg(j) = startSeq;
            if j < length(Y{v}.L)
                endSeq = Y{v}.seg(j+1)-1;
            else
                endSeq = Y{v}.seg(j+1);
                GT{v}.seg(j+1) = endSeq;
            end
            seq{v}(startSeq:endSeq,:) = X{v}(startSeq:endSeq,:);
            GT{v}.Lfr(startSeq:endSeq) = Y{v}.L(j);
        end
        if nSampGest > 0
            GTf = zeros(size(GT{v}.Lfr)); seqf = zeros(size(seq{v}));
            for k = 1:length(unique(GT{v}.Lfr))
                idx = find(GT{v}.Lfr == k);
                if isempty(idx)
                    if v == 1
                        devData = 'Training';
                    elseif v == 2
                        devData = 'Validation';
                    end
                    error('getTestSequences:NotEnoughSamples','There are no samples of gesture %d in the %s data',k,devData);
                end
                g = 1; j = 2;
                in = idx(1);
                while g <= nSampGest && j < length(idx)
                    if idx(j)-idx(j-1) > 1
                        fi = idx(j-1);
                        GTf(in:fi) = GT{v}.Lfr(in:fi);
                        seqf(in:fi,:) = seq{v}(in:fi,:);
                        in = idx(j);
                        g = g + 1;
                    end
                    j = j + 1;
                end            
            end
            idx = GTf ~= 0;
            GT{v}.Lfr = GT{v}.Lfr(idx);
            seq{v} = seq{v}(idx,:);
            SG = zeros(1,length(unique(Y{v}.L)));
            GT{v}.L = zeros(1,length(SG)*nSampGest);
            c = 1; k = 1; exceed = false;
            while ~all(SG == nSampGest) && ~exceed
                if c <= length(Y{v}.L)
                    if Y{v}.L(c) <= length(SG)
                        gesture = Y{v}.L(c);
                        if SG(gesture) < nSampGest
                            SG(gesture) = SG(gesture) + 1;
                            GT{v}.L(k) = Y{v}.L(c);
                            k = k + 1;
                        end
                    end
                    c = c + 1;
                else
                    exceed = true;
                end
            end
        else
            GT{v}.L = Y{v}.L;
        end
    else
        fps = 20;
        GTf = cell(1,length(l));
        seq{v} = cell(1,length(l));
        GT{v}.L = [];
        for i = 1:length(l)
            j = l(i);
            GTf{i} = zeros(1,fps*secs);
            seq{v}{i} = zeros(fps*secs,size(X{v},2));
            % Obtain continuous sequence of time fps*secs (approx)
            while Y{v}.seg(j) < Y{v}.seg(l(i))+fps*secs-1
                j = j + 1;
                if j > l(i)+1
                    startSeq = endSeq + 1;
                else
                    startSeq = 1;
                end    
                endSeq = startSeq + (Y{v}.seg(j)-1-Y{v}.seg(j-1));
                seq{v}{i}(startSeq:endSeq,:) = X{v}(Y{v}.seg(j-1):Y{v}.seg(j)-1,:);
                GTf{i}(startSeq:endSeq) = Y{v}.L(j-1);
                GT{v}.L = [GT{v}.L Y{v}.L(j-1)];
            end
        end
        GT{v}.Lfr = cell2mat(GTf); seq{v} = toMat(seq{v});
    end
end
    
%% Remove noise from training. We do not remove noise from validation
if ~noise
    idxDel = GT{1}.Lfr >= length(unique(Y{1}.L));
    GT{1}.Lfr(idxDel) = []; seq{1}(idxDel,:) = [];
    idxDel = GT{1}.L >= length(unique(Y{1}.L));
    GT{1}.L(idxDel) = [];
%     idxDel = GT{2}.Lfr >= length(unique(Y{2}.L));
%     GT{2}.Lfr(idxDel) = []; seq{2}(idxDel,:) = [];
%     idxDel = GT{2}.L >= length(unique(Y{2}.L));
%     GT{2}.L(idxDel) = [];
end
