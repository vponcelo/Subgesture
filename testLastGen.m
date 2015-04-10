function s = testLastGen(state,model,X,Y)

[Xdev,Ydev] = getDevSequences(X,Y,[],true,[],0);
Xval = Xdev{2}; Yval = Ydev{2};    
[~,bestPos] = min(state.Score);
I = state.Population(bestPos,:);
Iseg = round(I(2:end));
Iseg = Iseg(Iseg < inf);    % choose non-infinity segments
if model.probSeg > 0
    seg = round(reshape(Iseg,[2 size(Iseg,2)/2]));        
elseif model.probSeg == 0
    seg = round(reshape(I(2:end),[2 (size(I,2)-1)/2]));
end
if ~model.phmm.hmm                
    [~,s,~] = g(model,Xval,Yval);    % learn&optimize over validation
else
    if ~strcmp(model.phmm.clustType,'none')
        display('Discretizing validation sequence in Key Poses ...');
        Xval = discretizeSequence(model.Ctrain,Xval);
    end
    %% Discretize validation sequence
    display('Discretizing validation sequence in Subgestures ...');
    sw = model.sw;
    if sw > 0
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
        for sample = 1:length(Xval)
            KT = getUpdatedCosts(Xval{sample},model.SM);
            [~,Dval{sample}] = min(KT);
        end
    else
        KT = getUpdatedCosts(Xval,model.SM);
        [~,Dval] = min(KT);
    end
    [~,s] = evalswHMM(model, Dval, Yval, model.phmm.hmmTR, model.phmm.hmmE, model.phmm.model);
end