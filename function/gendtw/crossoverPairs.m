function xoverKids = crossoverPairs(parents,options,GenomeLength,~,~,thisPopulation,params)




% How many children to produce?
nKids = length(parents)/2;
% Extract information about linear constraints, if any
linCon = options.LinearConstr;
constr = ~isequal(linCon.type,'unconstrained');
% Allocate space for the kids
xoverKids = zeros(nKids,GenomeLength);

% To move through the parents twice as fast as thekids are
% being produced, a separate index for the parents is needed
index = 1;
% for each kid...
for i=1:nKids
    % get parents
    r1 = parents(index);
    index = index + 1;
    r2 = parents(index);
    index = index + 1;
    % Randomly select half of the genes from each parent
    % This loop may seem like brute force, but it is twice as fast as the
    % vectorized version, because it does no allocation.
    if strcmp(params.msmType,'none')
        ngenes = 0;  
    elseif strcmp(params.msmType,'fix')
        ngenes = 2;
    end
    for j = 1:GenomeLength-ngenes
        if j == 1
            k = 1;
        else
            k = j+1;
        end
        if(rand > 0.5)
            xoverKids(i,j:k) = thisPopulation(r1,j:k);
            
            if strcmp(params.msmType,'fix') && j == GenomeLength-ngenes
                % offspring of parents for the last 2 genes
                xoverKids(i,end-1:end) = thisPopulation(r1,end-1:end);
            end
        else
            xoverKids(i,j:k) = thisPopulation(r2,j:k);
            
            if strcmp(params.msmType,'fix') && j == GenomeLength-ngenes
                % offspring of parents for the last 2 genes
                xoverKids(i,end-1:end) = thisPopulation(r2,end-1:end);
            end
        end
    end    
    % Make sure that offspring are feasible w.r.t. linear constraints
    if constr
        feasible  = isTrialFeasible(xoverKids(i,:)',linCon.Aineq,linCon.bineq,linCon.Aeq, ...
            linCon.beq,linCon.lb,linCon.ub,sqrt(options.TolCon));
        if ~feasible % Kid is not feasible
            % Children are arithmetic mean of two parents (feasible w.r.t
            % linear constraints)
            alpha = rand;
            xoverKids(i,:) = alpha*thisPopulation(r1,:) + ...
                (1-alpha)*thisPopulation(r2,:);
        end
    end
end