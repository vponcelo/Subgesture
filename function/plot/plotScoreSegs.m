
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
title(sprintf('Mean scores throughout %d generations',currentGeneration));
% legend('Euclidean','Model');
xlabel('Number of segments');
ylabel('Mean scores');

global S; global Stest; global BESTIND;
try display(sprintf('Best k: %d\n',length(BESTIND(end).model.SM))); catch; end;

if length(S) > 1
    if state.Generation == 2 || S(end) > S(end-1) || Stest(end) >  Stest(end-1)
        global DATATYPE; global CACHE; global COORDS; global JOINTS; global NAT;
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
            params.Baseline,'_',params.mType,'_',num2str(currentGeneration),'gens','_',...
            num2str(length(JOINTS)),'joints',COORDS,'_','mod',num2str(NAT));
        try        
            save(strcat(filename,'.mat'),'S','Stest','CACHE','state','BESTIND','options','-v7.3');
            hgsave(gcf,filename,'-v7.3');
        catch e
            display(e.message);
        end
    end
end
if currentGeneration < params.generations
    bestind.model = []; bestind.state = []; BESTIND = [BESTIND bestind];
end