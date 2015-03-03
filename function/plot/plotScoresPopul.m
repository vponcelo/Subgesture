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
if strcmp(params.scoreMeasure,'overlap')
    title(sprintf('Mean population overlaps throughout %d generations',currentGeneration));
elseif strcmp(params.scoreMeasure,'levenshtein')
    title(sprintf('Mean levenshtein distances throughout %d generations',currentGeneration));
end
% legend('Euclidean','Model');
xlabel('Population');
if strcmp(params.scoreMeasure,'overlap')
    ylabel('Mean overlaps');
elseif strcmp(params.scoreMeasure,'levenshtein')
    ylabel('Mean levenshtein distances');
end
