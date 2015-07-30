function s = testLastGen(state,model,X,Y,X_l)

if ~iscell(X)
    Y.Lfr = [];
    for i = 1:length(Y.seg)-1
        Y.Lfr = [Y.Lfr repmat(Y.L(i),1,Y.seg(i+1)-Y.seg(i))];
        if i == length(Y.seg)-1, Y.Lfr = [Y.Lfr Y.Lfr(end)]; end
    end
else
%     [Xdev,Ydev] = getDevSequences(X,Y,[],true,[],0);
%     Xval = Xdev{2}; Yval = Ydev{2};
end
if model.classification
    Xval = X_l;
else
    Xval = X;
end
Yval = Y;

%% Decode Best Individual
[~,bestPos] = min(state.Score);
I = state.Population(bestPos,:);
Iseg = round(I(2:end));
Iseg = Iseg(Iseg < inf);    % choose non-infinity segments
if model.probSeg > 0
    seg = round(reshape(Iseg,[2 size(Iseg,2)/2]));        
elseif model.probSeg == 0
    seg = round(reshape(I(2:end),[2 (size(I,2)-1)/2]));
end

%% Evaluation
if ~model.phmm.hmm
    if model.darwin
        display('Obtaining Validation Sequence from Subgestures U ...')
        [~,model.KT] = computeSGFromU(model.Us,Xval);
    end
    [~,s,~] = g(model,Xval,Yval);    % learn&optimize over validation
else
    if ~strcmp(model.phmm.clustType,'none')
        display('Discretizing validation sequence in Key Poses ...');
        Xval = discretizeSequence(model.Ctrain,Xval);
    end
    %% Discretize validation sequence
    display('Discretizing validation sequence in Subgestures ...');
    sw = model.sw;
    if sw > 0 && ~model.classification
        r = inf;
        while any(r > length(Xval)-sw)
            r=randperm(round(length(Xval)),1);
        end
        sseg=r(1):min(r(1)+sw,length(Yval.Lfr));
        Xval=Xval(sseg,:);
        Yval.Lfr=Yval.Lfr(sseg);
    end
    Dval = cell(1,length(Xval));
    if iscell(Xval)
        if model.darwin
            Dval = computeSGFromU(model.Us,Xval);
        else
            Dval = cell(length(Xval));
            for l = 1:length(Xval)
                Dval{l} = cell(length(Xval{l}));
                for sample = 1:length(Xval{l})
                    KT = getUpdatedCosts(Xval{l}{sample},model.SM);
                    [~,Dval{l}{sample}] = min(KT);
                end
            end
        end
    else
        if model.darwin
            Dval = computeSGFromU(model.Us,Xval);
        else
            KT = getUpdatedCosts(Xval,model.SM);
            [~,Dval] = min(KT);
        end
    end
    [~,s] = evalswHMM(model, Dval, Yval);
end