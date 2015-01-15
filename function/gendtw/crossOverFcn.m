function xoverKids = crossOverFcn(parents,options,nvars,FitnessFcn,unused,thisPopulation,params)

%% Crossover
if rand() >= params.thMutCross
    xoverKids = crossoverPairs(parents,options,nvars,FitnessFcn,unused,thisPopulation);
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
    Peven = xoverKids(i,2:2:end);
    Podd = xoverKids(i,3:2:end);
    Psort = reshape(sortrows([Peven; Podd]')',1,size(xoverKids,2)-1);
    xoverKids(i,2:end) = Psort;
    
    cutPoint = find(xoverKids(i,:) == inf,1);
    if mod(cutPoint,2) > 0
        cutPoint = cutPoint - 1;
    end
    xoverKids(i,cutPoint:end) = inf;    
end
