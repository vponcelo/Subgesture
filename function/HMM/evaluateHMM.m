function [prob] = evaluateHMM(seqTest, TRANS, EMIS, modelpmtk)
%
% Donada una sequencia de test 'seqTest' i un model apres (matrius 'TRANS',
% 'EMIS'), retorna una mesura de probabilitat que 'seqTest' s'hagi generat
% amb el model apres
%

% S'estimen els estats versemblants:
try
    if isempty(modelpmtk)
        likelystates = hmmviterbi(seqTest, TRANS, EMIS);
        TR = TRANS;
    else
        likelystates = hmmMap(modelpmtk, seqTest);
%         likelyStates = hmmviterbi(seqTest, modelpmtk.transPrior, modelpmtk.emission.T);   % Equivalent
%         [mViterbiErr,mMaxMargErr,emViterbiErr, ...
%             emMaxMargErr] = decodeHMMStates(modelpmtk, likelystates, hidden, seqTest);
        TR = modelpmtk.A;        
    end

    % Es calcula la mesura de versemblança
    prob=1;
    len=length(likelystates)-1;
    for i=1:len
        e1=likelystates(i);
        e2=likelystates(i+1);
        prob=prob*TR(e1,e2);
    end
catch error
    % No path in the sequence
    warning(error.identifier,error.message);
%     if strcmp(error.identifier,'stats:hmmviterbi:ZeroTransitionProbability'),
%         prob=0;
%     end
    prob=0;
end

