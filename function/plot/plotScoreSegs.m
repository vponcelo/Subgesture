function state = plotScoreSegs(options,state,~,params)
% plot mean scores

if state.Generation < 1
    return;
end

if length(state.Best) ~= state.Generation
    error('There must be as scores as generations');
end

x = sum(state.Population(:,2:end)' < inf)/2;

if strcmp(params.scoreMeasure,'overlap');
    s = -state.Score;    
elseif strcmp(params.scoreMeasure,'levenshtein');
    s = state.Best(end);    
end
s(s < 0) = 0;

plot(x,s,'o');
if strcmp(params.scoreMeasure,'overlap')
    title(sprintf('Mean overlaps throughout %d generations',params.generations-options.Generations+state.Generation));
elseif strcmp(params.scoreMeasure,'levenshtein')
    title(sprintf('Mean levenshtein distances throughout %d generations',params.generations-options.Generations+state.Generation));
end
% legend('Euclidean','Model');
xlabel('Number of segments');
if strcmp(params.scoreMeasure,'overlap')
    ylabel('Mean overlap');
elseif strcmp(params.scoreMeasure,'levenshtein')
    ylabel('Mean levenshtein distances');
end
