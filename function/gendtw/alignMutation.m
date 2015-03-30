function P = alignMutation(parents,state,options,~,thisPopulation,params,seg0)
% Output:
%   P: Mutation Children
% Input:
%   parents: Population parents
%   GenomeLength: length of the genome (nvars)
%   thisPopulation: current Population
%   seg: Initial derivative segmenetation

global BASELINE;

if (params.shrink > 1) || (params.shrink < 0)
    warning(message('globaloptim:mutationgaussian:shrinkFactor'));
end

scale = params.scale - params.shrink * params.scale * state.Generation/options.Generations;

range = options.PopInitRange;
lower = range(1,:);
upper = range(2,:);
scale = scale * (upper - lower);
    
if ~strcmp(params.Baseline,BASELINE{1})
    seg0ini = seg0(1:2:end);
    seg0fi = seg0(2:2:end);
end
  
P = thisPopulation(parents,:);

for i=1:length(parents)
    parent = P(i,1:params.N*2+1);   % In the case of fixed MSM, the last 2 genes don't mutate here
    for j = 1:length(parent)
        if j == 1
            P(i,j) = parent(j)  + scale(j) .* randn();
%             P(i,j) = parent(j);
        elseif mod(j,2) == 0
            if P(i,j) < inf
                if strcmp(params.Baseline,BASELINE{2})
                    % generate derivate population
                    [~,pos] = min(abs(seg0ini-parent(j)));
                    P(i,j) = seg0ini(pos(1));
                end
            end
        else
            if P(i,j) < inf
                if strcmp(params.Baseline,BASELINE{2})
                    % generate derivate population 
                    [~,pos] = min(abs(seg0fi-(parent(j-1)+parent(j)-1)));
                    if seg0fi(pos(1)) < parent(j-1)
                        if pos(1)+1 > length(seg0fi);
                            [~,posIn] = min(abs(seg0ini-parent(j-1)));
                            fi = seg0ini(posIn);
                        else
                            fi = seg0fi(pos(1)+1);
                        end
                    else
                        fi = seg0fi(pos(1));
                    end
                    P(i,j) = fi-parent(j-1)+1;
                end
            end
        end
    end
    if strcmp(params.mType,'modelSM1') || strcmp(params.mType,'allSM1')
        P(i,end-1) = P(i,end-1)  + scale(end-1) .* randn();
        P(i,end) = P(i,end)  + scale(end) .* randn();
    end
end
