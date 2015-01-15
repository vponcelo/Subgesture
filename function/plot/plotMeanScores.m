function state = plotMeanScores(options,state,~,params,S_eu)
% plot mean scores
% S_eu: Mean score directly obtained from DTW with euclidean distance

if state.Generation < 1
    warning('Current generation is lesser than 1');
    return;
end

if length(state.Best) ~= state.Generation
    error('There must be as scores as generations');
end

fprintf('Generation %d\n\n',params.generations-options.Generations+state.Generation);

global S;

length(S)

x = 1:1:params.generations-options.Generations+state.Generation;
length(x)
S_eu = S_eu*ones(1,length(x));

if strcmp(params.scoreMeasure,'overlap');
    s = -state.Best(end);
    if state.Generation > 1
        if s < -state.Best(end-1)
            s = -state.Best(end-1); 
        end
    end
elseif strcmp(params.scoreMeasure,'levenshtein');
    s = state.Best(end);
    if state.Generation > 1
        if s > state.Best(end-1)
            s = state.Best(end-1); 
        end
    end
end
if s < 0
    s = 0;
end

if length(S) >= options.Generations
    warning('Score length is greater or equal to the total number of generations');
    return;
end

S = [S s];

length(S)

if length(S) ~= length(x)
    warning('Score length is must be equal to the x axis in order to plot');
	return;
end

plot(x,S_eu,'k');
hold on
plot(x,S,'b');
hold off
if strcmp(params.scoreMeasure,'overlap')
    title(sprintf('Mean overlaps throughout %d generations',params.generations-options.Generations+state.Generation));
elseif strcmp(params.scoreMeasure,'levenshtein')
    title(sprintf('Mean levenshtein distances throughout %d generations',params.generations-options.Generations+state.Generation));
end
% legend('Euclidean','Model');
xlabel('Generation');
if strcmp(params.scoreMeasure,'overlap')
    ylabel('Mean overlaps');
elseif strcmp(params.scoreMeasure,'levenshtein')
    ylabel('Mean levenshtein distances');
end

if state.Generation > 0 && mod(state.Generation,1) == 0
    global DATATYPE;
    global PERCENTDATA;
    global COORDS;
    global JOINTS;
    global NAT;
%     if strcmp(params.scoreMeasure,'overlap')
        filename = strcat('results/',DATATYPE,'/validation/Exp3/',params.Baseline,'Results',num2str(params.generations-options.Generations+state.Generation),'_',num2str(length(JOINTS)),COORDS,num2str(PERCENTDATA),'%_',num2str(NAT));
        try
            save(strcat(filename,'.mat'),'S','CACHE','state','options','MODEL','-v7.3');
            hgsave(gcf,filename);
        catch e
            %display(e.Message);
        end
%     elseif strcmp(params.scoreMeasure,'levenshtein')
%         filename = strcat('results/',DATATYPE,'/validation/Exp3/',params.Baseline,'Results',num2str(state.Generation),'_',num2str(PERCENTDATA),'%');
%         save(strcat(filename,'.mat'),'S','CACHE','state','options','MODEL','-v7.3');
%         hgsave(gcf,filename);
%     end
end