function [trueModelViterbiErr,trueModelMaxMargErr,emModelViterbiErr,emModelMaxMargErr] = decodeHMMStates(model, path, hidden, observed)
%% Decode using true model
try
    decodedFromTrueViterbi = bestPermutation(path, hidden);
    trueModelViterbiErr = mean(decodedFromTrueViterbi ~= hidden);
catch
end
    decodedFromTrueMaxMarg = maxidx(hmmInferNodes(model, observed), [], 1);
try
    decodedFromTrueMaxMarg = bestPermutation(decodedFromTrueMaxMarg, hidden);
    trueModelMaxMargErr = mean(decodedFromTrueMaxMarg ~= hidden);
catch
end

%% Decode using the EM model
try
    modelEM = hmmFit(observed, params.phmm.states, 'discrete', ...
            'convTol', 1e-5, 'nRandomRestarts', 2, 'verbose', false);
    decodedFromEMviterbi = hmmMap(modelEM, observed);
    decodedFromEMviterbi = bestPermutation(decodedFromEMviterbi, hidden);
    emModelViterbiErr = mean(decodedFromEMviterbi ~= hidden);
    decodedFromEMmaxMarg = maxidx(hmmInferNodes(modelEM, observed), [], 1);
    decodedFromEMmaxMarg = bestPermutation(decodedFromEMmaxMarg, hidden);
    emModelMaxMargErr = mean(decodedFromEMmaxMarg ~= hidden);
catch
end