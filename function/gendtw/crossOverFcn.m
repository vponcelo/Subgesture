function xoverKids = crossOverFcn(parents,options,nvars,FitnessFcn,unused,thisPopulation,params,X)

%% Crossover
if rand() >= params.probCross
    xoverKids = crossoverPairs(parents,options,nvars,FitnessFcn,unused,thisPopulation,params);
else
    xoverKids = crossoverscattered(parents,options,nvars,FitnessFcn,unused,thisPopulation);
end

%% Sort xoverKids and add invalid segments
for i = 1:size(xoverKids,1)
    %% Old implementation
%     [xoverKids(i,2:2:end),is] = sort(xoverKids(i,2:2:end),'ascend');
%     lengths = xoverKids(i,3:2:end);
%     for j = 1:length(is)
%         xoverKids(i,2*j+1) = lengths(is(j));
%         if j > 1
%             if xoverKids(i,2*j) == xoverKids(i,2*(j-1))
%                 if xoverKids(i,2*j+1) < xoverKids(i,2*(j-1)+1)
%                     next = xoverKids(i,2*(j-1)+1);
%                     xoverKids(i,2*(j-1)+1) = xoverKids(i,2*j+1);
%                     xoverKids(i,2*j+1) = next;
%                 end
%             end
%         end
%     end

    %% this run faster
    % Sort P
    Peven = xoverKids(i,2:2:params.N*2);
    Podd = xoverKids(i,3:2:params.N*2+1);
    Psort = reshape(sortrows([Peven; Podd]')',1,params.N*2);
    xoverKids(i,2:params.N*2+1) = Psort;
    
    cutPoint = find(xoverKids(i,:) == inf,1);   % find the first infinity value
    if mod(cutPoint,2) > 0
        cutPoint = cutPoint - 1;
    end
    xoverKids(i,cutPoint:params.N*2+1) = inf;    
end

xoverKids = repairFcn(xoverKids, params, length(X), nvars);
