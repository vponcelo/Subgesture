function state = plotScoreSegs(options,state,~,params)
% plot mean scores

if state.Generation < 1
    return;
end

x = sum(state.Population(:,2:end)' < inf)/2;

if strcmp(params.scoreMeasure,'overlap');
    s = -state.Score;    
elseif strcmp(params.scoreMeasure,'levenshtein');
    s = state.Best(end);    
end
s(s < 0) = 0;

currentGeneration = params.generations-options.Generations+state.Generation;

plot(x,s,'o');
if strcmp(params.scoreMeasure,'overlap')
    title(sprintf('Mean overlaps throughout %d generations',currentGeneration));
elseif strcmp(params.scoreMeasure,'levenshtein')
    title(sprintf('Mean levenshtein distances throughout %d generations',currentGeneration));
end
% legend('Euclidean','Model');
xlabel('Number of segments');
if strcmp(params.scoreMeasure,'overlap')
    ylabel('Mean overlap');
elseif strcmp(params.scoreMeasure,'levenshtein')
    ylabel('Mean levenshtein distances');
end

global S;

if length(S) > 1
    if state.Generation > 0 && S(end) > S(end-1) %&& mod(currentGeneration,1) == 0
        global DATATYPE;
        global CACHE;
        global MODEL;
        global PERCENTDATA;
        global COORDS;
        global JOINTS;
        global NAT;
        if ~exist(strcat('results/',DATATYPE,'/validation/Exp3/gen',num2str(params.generations),'popul',num2str(options.PopulationSize)),'dir')
            mkdir(strcat('results/',DATATYPE,'/validation/Exp3/gen',num2str(params.generations),'popul',num2str(options.PopulationSize)));
        end
        filenames = dir(strcat('results/',DATATYPE,'/validation/Exp3/gen',num2str(params.generations),'popul',num2str(options.PopulationSize),'/'));
        if ~isempty(filenames)
            filename = strcat('results/',DATATYPE,'/validation/Exp3/gen',num2str(params.generations),'popul',num2str(options.PopulationSize),'/',filenames(end).name(1:end-4));
            if exist(strcat(filename,'.mat'),'file')
                delete(strcat(filename,'.mat'));
            end
            if exist(strcat(filename,'.fig'),'file')
                delete(strcat(filename,'.fig'));
            end
        end
        set(gcf, 'Position', [0 0 1920 1200]);
        filename = strcat('results/',DATATYPE,'/validation/Exp3/gen',num2str(params.generations),'popul',num2str(options.PopulationSize),'/',...
            params.Baseline,'_',params.msmType,'_',num2str(currentGeneration),'gens','_',...
            num2str(length(JOINTS)),'joints',COORDS,'_','mod',num2str(NAT));
        try        
            save(strcat(filename,'.mat'),'S','CACHE','state','options','MODEL','-v7.3');
            hgsave(gcf,filename,'-v7.3');
        catch e
            display(e.message);
        end
    end
end