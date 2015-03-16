function [Z,C] = kmeansDTW_v0(X,k,maxIter,nConverg)
% K-means DTW algorithm that obtains k subsequences from X
%
% Input:
%   X: Sequence Data divided into cells
%   k: number of initial clusters (data subsets)
%   maxIter: number of maximum iterations (default:20)
%   nConverg: number of repetitions to converge (default:2)
% Output:
%   C: array with the cluster indices for which each sequence of X belongs to
%   Z: cell with the k centroide sequences

if ~iscell(X) && length(X) > k
    error('X must be a cell with length greater than k containing a vector at each position');
end
if k <= 1
    error('k must be a natural number greater than 1 (minimum of 2 data subsets)');
end

if ~exist('maxIter','var')
    maxIter = 20;
elseif maxIter < 1
    error('number of repetitions to converge must be a natural number greater than 0.');
end

if ~exist('nConverg','var')
    nConverg = 2;
elseif nConverg < 1
    error('number of iterations must be a natural number greater than 0.');
end


%% randomly obtain initial (and different) centroids
Z = cell(1,k);
Zi = zeros(1,k);
Zi(1) = randi(length(X));
for i = 2:k
    Zi(i) = randi(length(X));
    for j = i-1:-1:1
        while Zi(i) == Zi(j)
            Zi(i) = randi(length(X));
        end
    end
end    

%% Initialize centroid sequences and iterate until convergence or maximum Iterations reached
dc = inf*ones(k,length(X));
dc2 = inf*ones(k,length(X));
count = nConverg;
it = 0;
while it < maxIter && count > 0
    for j = 1:k 
        for i = 1:length(X)
            if i ~= Zi(j)
                W=dtwc(X{i},X{Zi(j)},1);
                dc(j,i) = W(end,end);
            else
                dc(j,i) = 0;
            end
        end
    end 
    text = '';
    convergs = false(1,k);
    for i = 1:k 
        dc(i,Zi(i)) = inf;
        if min(dc(i,:)) <= min(dc2(i,:))
            [cost,mincPos] = min(dc(i,:));
            for j = 1:k
                if j ~= i                    
                    if mincPos == Zi(j)
                        dc(i,mincPos) = inf;
                        [~,mincPos] = min(dc(i,:));
                    end                    
                end
            end
            if dc(i,mincPos) == min(dc2(i,:))
                convergs(i) = true;
                if all(convergs)
                    count = count - 1;
                end
            else
                convergs(i) = false;
                count = nConverg;
            end
        else
            [~,mincPos] = min(dc2(i,:));
            dc(i,:) = dc2(i,:);
            convergs(i) = false;            
            count = nConverg;
        end            
        Zi(i) = mincPos;
        Z{i} = X{Zi(i)};
%         text = strcat(text,num2str(dc(i,c(i))),'\t');        
    end         
%     fprintf(sprintf('%s\n',text));
    dc2 = dc;
    it = it + 1;  
end

%% Save final subsets
C = zeros(1,length(X));
for i = 1:length(X)
    [~,kf] = min(dc(:,i));
    C(i) = kf;
end
