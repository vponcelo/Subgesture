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

display(sprintf('Generation %d\n\n',params.generations-options.Generations+state.Generation));

x = 1:1:params.generations-options.Generations+state.Generation;
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
% if s < 0
%     s = 0;
% end

global S;

if length(S) >= options.Generations
    warning('Score length is greater or equal to the total number of generations');
    return;
end

try 
    s_end = S(end);
    if s < s_end
        s = s_end;
    end
catch
end

S = [S s];

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
    global CACHE;
    global MODEL;
    global PERCENTDATA;
    global COORDS;
    global JOINTS;
    global NAT;
%     if strcmp(params.scoreMeasure,'overlap')
        if ~exist(strcat('results/',DATATYPE,'/validation/Exp3/gen',num2str(options.Generations),'popul',num2str(options.PopulationSize)),'dir')
            mkdir(strcat('results/',DATATYPE,'/validation/Exp3/gen',num2str(options.Generations),'popul',num2str(options.PopulationSize)));
        end
        filename = strcat('results/',DATATYPE,'/validation/Exp3/gen',num2str(options.Generations),'popul',num2str(options.PopulationSize),'/',...
            params.Baseline,'_',params.msmType,'_',num2str(params.generations-options.Generations+state.Generation),'gens','_',...
            num2str(length(JOINTS)),'joints',COORDS,'_','mod',num2str(NAT));
        try
            save(strcat(filename,'.mat'),'S','CACHE','state','options','MODEL','-v7.3');
            set(gcf, 'Position', [0 0 1920 1200]);
            hgsave(gcf,filename,'-v7.3');
            %saveas(gcf,strcat(filename,'Resized'),'png');            
            filename = strcat('results/',DATATYPE,'/validation/Exp3/gen',num2str(options.Generations),'popul',num2str(options.PopulationSize),'/',...
                params.Baseline,'_',params.msmType,'_',num2str(params.generations-options.Generations+state.Generation-1),'gens','_',...
                num2str(length(JOINTS)),'joints',COORDS,'_','mod',num2str(NAT));
            if exist(strcat(filename,'.mat'),'file')
                delete(strcat(filename,'.mat'));
            end
            if exist(strcat(filename,'.fig'),'file')
                delete(strcat(filename,'.fig'));
            end
        catch e
            display(e.message);
        end
%     elseif strcmp(params.scoreMeasure,'levenshtein')
%         filename = strcat('results/',DATATYPE,'/validation/Exp3/',params.Baseline,'Results',num2str(state.Generation),'_',num2str(PERCENTDATA),'%');
%         save(strcat(filename,'.mat'),'S','CACHE','state','options','MODEL','-v7.3');
%         hgsave(gcf,filename);
%     end
end