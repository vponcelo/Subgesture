function state = plotScoresPopul(options,state,~,params)
% plot mean scores

if state.Generation < 1
    return;
end

if length(state.Best) ~= state.Generation
    error('There must be as scores as generations');
end

x = 1:1:params.population;

if strcmp(params.scoreMeasure,'overlap');
    s = -state.Score;    
elseif strcmp(params.scoreMeasure,'levenshtein');
    s = state.Score;    
end
s(s < 0) = 0;

currentGeneration = params.generations-options.Generations+state.Generation;

bar(x,s);
title(sprintf('Mean population scores throughout %d generations',currentGeneration));
% legend('Euclidean','Model');
xlabel('Population');
ylabel('Mean scores');
