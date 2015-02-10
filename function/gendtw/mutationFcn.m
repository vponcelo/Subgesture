function mutationChildren = mutationFcn(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation,params,seg0,X)

%% Get parents an add/delete segments given a probability
P = thisPopulation(parents,:);
for i = 1:length(parents)
    if params.probSeg > 0
        if strcmp(params.msmType,'fix')            
            randSegs = rand(1,(nvars-3)/2);
        else
            randSegs = rand(1,(nvars-1)/2);
        end
        
        %% this run faster        
%         addDel = randSegs <= params.probSeg;
%         idxMod = false(1,size(P,2));
%         for j = 2:length(idxMod)
%             idxMod(j) = addDel(round((j-1)/2));
%         end
%         idxInf = mod((P(i,idxMod) == inf),2) == 0;        
%         P(i,idxInf) = randi([1 length(X)-params.nmax+1],1,sum(idxInf));
%         idxInf = mod((P(i,idxMod) == inf),2) == 1;
%         P(i,idxInf) = randi([params.nmin params.nmax],1,sum(idxInf));

        %% Old implementation        
        for j = 2:params.N*2+1
            % mutate (add/delete) parents' segments if chosen
            if randSegs(round((j-1)/2)) <= params.probSeg
                if mod(j,2) == 0
                    if P(i,j) == inf
                        % generate random population
                        P(i,j) = randi([1 length(X)-params.nmax+1]);
                    else
                        P(i,j) = inf;
                    end
                else
                    if P(i,j) == inf
                         % generate random population
                        P(i,j) = randi([params.nmin params.nmax]);
                    else
                        P(i,j) = inf;
                    end
                end
            end
        end
        
        %% this run faster
        % Sort P
        Peven = P(i,2:2:params.N*2);
        Podd = P(i,3:2:params.N*2+1);
        Psort = reshape(sortrows([Peven; Podd]')',1,params.N*2);
        P(i,2:params.N*2+1) = Psort;
        
        %% Old implementation
%         [P(i,2:2:end),is] = sort(P(i,2:2:end),'ascend');
%         lengths = P(i,3:2:end);
%         for j = 1:length(is)
%             P(i,2*j+1) = lengths(is(j));
%             if j > 1
%                 if P(i,2*j) == P(i,2*(j-1))
%                     if P(i,2*j+1) < P(i,2*(j-1)+1)
%                         next = P(i,2*(j-1)+1);
%                         P(i,2*(j-1)+1) = P(i,2*j+1);
%                         P(i,2*j+1) = next;
%                     end
%                 end
%             end
%         end
%         nsegs = sum(transpose(P(i,:) < inf))+1;
%         if mod(nsegs,2) > 0
%             nsegs = nsegs - 1;
%         end
%         if P(i,1) > nsegs(i)
%             P(i,1) = randi([params.k0 nsegs(i)]);            
%         end
    end
end
thisPopulation(parents,:) = P;

%% Mutation type given the probability
if rand() >= params.thMutCross    
    mutationChildren = alignMutation(parents,state,options,nvars,thisPopulation,params,seg0);
else
    mutationChildren = mutationgaussian(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation,params.scale,params.shrink);
end

%% Sort mutationChildren
for i = 1:length(parents)
    %% this run faster
     % Sort mutationChildren
    Peven = mutationChildren(i,2:2:params.N*2);
    Podd = mutationChildren(i,3:2:params.N*2+1);
    Psort = reshape(sortrows([Peven; Podd]')',1,params.N*2);
    mutationChildren(i,2:params.N*2+1) = Psort;
    
    %% Old implementation
%     mutationChildren = round(mutationChildren);
%     [mutationChildren(i,2:2:end),is] = sort(mutationChildren(i,2:2:end),'ascend');
%     lengths = mutationChildren(i,3:2:end);
%     for j = 1:length(is)
%         mutationChildren(i,2*j+1) = lengths(is(j));
%         if j > 1
%             if mutationChildren(i,2*j) == mutationChildren(i,2*(j-1))
%                 if mutationChildren(i,2*j+1) < mutationChildren(i,2*(j-1)+1)
%                     next = mutationChildren(i,2*(j-1)+1);
%                     mutationChildren(i,2*(j-1)+1) = mutationChildren(i,2*j+1);
%                     mutationChildren(i,2*j+1) = next;
%                 end
%             end
%         end
%     end    
end
