function [Z,C] = kmeansDTW(X,k,version,distance,labelsT,maxIter)
% K-means DTW algorithm that obtains k subsequences from X
%
% Input:
%   X: Sequence Data divided into cells
%   k: number of initial clusters (data subsets)
%   version: string with the version of the kmeans DTW algorithm
%       % Possibilities: 'vM_I_E' M: Maximization type; I: Initialization
%       type; E: Expectation type (0 if is not included)
%           'v1_0': Path sum of minimum costs (1); random init (0)
%           'v1_1': Path sum of minimum costs (1); median sequences of each sorted k groups (1)
%           'v1_2': Path sum of minimum costs (1); k median sequences of sorted sequences (2)
%           'v1_3': Path sum of minimum costs (1); k median sequences of labeled gesture sequences (3) 
%           'v2_0_0': mean sequence alignment (2); random init (0);
%           Centroid CM (0) - NOT WORKING AT ALL: It remains to deal with empty centroids
%           'v2_0_1': mean sequence alignment (2); random init (0); Centroid closest sequence to CM (1)
%           'v2_1_0': mean sequence alignment (2); median sequences of each sorted k groups (1); Centroid CM (0)
%           'v2_1_1': mean sequence alignment (2); median sequences of each sorted k groups (1); Centroid closest sequence to CM (1)
%           'v2_2_0': mean sequence alignment (2); k median sequences of sorted sequences (2); Centroid CM (0)
%           'v2_2_1': mean sequence alignment (2); k median sequences of sorted sequences (2); Centroid closest sequence to CM (1)
%           'v2_3_0': mean sequence alignment (2); k median sequences of labeled gesture sequences (3); Centroid CM (0)
%           'v2_3_1': mean sequence alignment (2); k median sequences of labeled gesture sequences (3); Centroid closest sequence to CM (1)
%   distance: string with the type of distance
%           'dtwCost': DTW distance
%           'euclidean': euclidean distance
%   labelsT: labels for training data
%   maxIter: number of maximum iterations (default:50)
% Output:
%   C: array with the cluster indices for which each sequence of X belongs to
%   Z: cell with the k centroide sequences

if ~iscell(X) && length(X) > k
    error('X must be a cell with length greater than k containing a vector at each position');
end
if k <= 1
    error('k must be a natural number greater than 1 (minimum of 2 data clusters)');
end

if ~exist('maxIter','var')
    maxIter = 20;
elseif maxIter < 1
    error('number of repetitions to converge must be a natural number greater than 0.');
end

it = 0;
Z = cell(1,k);
convergs = false(1,k);
costsZ = -1*ones(1,k);
costs2Z = -1*ones(1,k);
%% Validation: Save expectations at each iteration
% V = -1*ones(k,maxIter);
% CG = zeros(length(X),maxIter);

%% Initialize centroid sequences 
if strcmp(version(4),'0')
    Zi = zeros(1,k);    
    for j = 1:k 
        aux = Zi;    
        while ismember(Zi(j),aux)
            Zi(j) = randi(length(X));
        end
        Z{j} = X{Zi(j)};    % Initialization
    end
elseif strcmp(version(4),'1')
    % obtain initial median sequences of each sorted k groups (1)
    Xsort = cell(1,length(X));
    slengths = zeros(1,length(X));
    for i = 1:length(X)
        slengths(i) = size(X{i},1);
    end
    for i = 1:length(X)
        [~,pos] = min(slengths);
        Xsort{i} = X{pos};
        slengths(pos) = max(slengths)+1;
    end
    nseqsGroup = round(length(X)/k);
    aux = Xsort;
    for j = 1:k
        slengths = zeros(1,nseqsGroup);
        i = 1;
        while i <= nseqsGroup && ~isempty(aux)
            slengths(i) = size(aux{1},1);
            aux(1) = [];
            i = i + 1;
        end
        idx = slengths == 0;
        slengths(idx) = [];
        if mod(length(slengths),2) == 0
            meanLength = mean(slengths);
            dists2mean = abs(slengths-meanLength);
            [~,idxm] = min(dists2mean);
        else
            medianLength = median(slengths);
            idxm = find(slengths == medianLength);
        end   
        if length(unique(idxm)) > 1
            idxm = idxm(round(length(idxm)/2));
        end              
        Z{j} = Xsort{(j-1)*nseqsGroup+idxm};     % Initialization
    end
elseif strcmp(version(4),'2')
    % obtain initial k median sequences of sorted sequences
    slengths = zeros(1,length(X));
    for i = 1:length(X)
        slengths(i) = size(X{i},1);
    end
    aux = slengths;
    idx = true(1,length(aux));
    for j = 1:k
        if mod(length(aux(idx)),2) == 0
            meanLength = mean(aux(idx));
            dists2mean = abs(aux-meanLength);
            dists2mean(~idx) = inf;
            [~,idxm] = min(dists2mean);                
        else
            medianLength = median(aux(idx));
            idxm = find(aux == medianLength);
            idxm(ismember(idxm,find(~idx))) = [];
        end
        if length(unique(idxm)) > 1
            idxm = idxm(round(length(idxm)/2));
        end
        if ~isscalar(idxm)  % should not happen
            error('More than one median sequence detected. Check indices.');
        end     
        idx(idxm) = false;
        Z{j} = X{idxm};     % Initialization
    end
    clear aux;
elseif strcmp(version(4),'3')
    % obtain initial k median sequences of labeled gesture sequences 
    for i = 1:k
        X_k = X(labelsT == i);
        Xsort = cell(1,length(X_k));
        slengths = zeros(1,length(Xsort));
        for j= 1:length(Xsort)
            slengths(j) = size(X_k{j},1);
        end
        aux = slengths;
        for j = 1:length(Xsort)
            [~,pos] = min(aux);
            Xsort{j} = X_k{pos};
            aux(pos) = max(aux)+1;
        end
        clear aux;            
        if mod(length(slengths),2) == 0
            meanLength = mean(slengths);
            dists2mean = abs(slengths-meanLength);
            [~,idxm] = min(dists2mean);
        else
            medianLength = median(slengths);
            idxm = find(slengths == medianLength);
        end   
        if length(unique(idxm)) > 1
            idxm = idxm(round(length(idxm)/2));
        end
        if ~isscalar(idxm)  % should not happen
            error('More than one median sequence detected. Check indices.');
        end
        Z{i} = Xsort{idxm};     % Initialization
    end
end

%% Main loop: iterate until convergence or maximum iterations reached
while it < maxIter && ~all(convergs) 
    Z2 = Z;
    %% Maximization (update the centroids)
    if it > 0
%         text = sprintf('Iter %d centroid distances: ',it);
        for j = 1:k 
            idx = find(C == j);
            if isempty(idx)
                Z{j} = [];
                costsZ(j) = -1;
                if isequal(Z{j},Z2{j})
                    convergs(j) = true;
                end
                costs2Z(j) = costsZ(j);                
                continue;
            end
            if strcmp(version(2),'1')
                %% Assign the minimum cost sums of DTW matrices between all current cluster sequence combinations
                dcf = -1*ones(1,length(idx));
                for i = 1:length(idx)
                    dc = zeros(1,length(idx)-1);
                    for l = 1:length(idx)
                        if i ~= l
                            if strcmp(distance,'dtwCost')
                                W=dtwc(X{idx(i)},X{idx(l)},1);
                                dc(l) = W(end,end);
                            elseif strcmp(distance,'euclidean')   % only possible if sequence lenghts are the same
                                dc(l) = sum(diag(pdist2(X{idx(i)},X{idx(l)})))/length(X{idx(i)});
                            end
                        end
                    end
                    dcf(i) = sum(dc);
                end
                [~,iz] = min(dcf);               
                Z{j} = X{idx(iz)};
            elseif strcmp(version(2),'2')
                %% Obtain the mean aligned sequence
                % Select median sequence
                slengths = zeros(1,length(X(idx)));        
                for i = 1:length(slengths)
                    slengths(i) = size(X{idx(i)},1);
                end        
                if mod(length(slengths),2) == 0
                    meanLength = mean(slengths);
                    dists2mean = abs(slengths-meanLength);
                    [~,idxm] = min(dists2mean);
                else
                    medianLength = median(slengths);
                    idxm = find(slengths == medianLength);
                end   
                if length(unique(idxm)) > 1
                    idxm = idxm(round(length(idxm)/2));
                end
                ptr = X{idx(idxm)};

                % Compute the mean alignment respect to the median sequence
                alig_seqs = zeros(length(idx),size(ptr,1),size(ptr,2));
                for i = 1:length(idx)
                    if i ~= idx
                        if strcmp(distance,'dtwCost')
                            W=dtwc(X{idx(i)},ptr,1);
                            [~,~,A] = aligngesture(X{idx(i)},W);
                            if size(A,1) ~= size(alig_seqs,2)
                                error('All elements of the cost matrix W are 0');
                            end
                            alig_seqs(i,:,:) = A;
                        elseif strcmp(distance,'euclidean')
                            alig_seqs(i,:,:) = X{idx(i)};
                        end
                    else
                        alig_seqs(i,:,:)=ptr;
                    end
                end 
                if size(alig_seqs,1) > 1
                    CM = reshape(mean(alig_seqs),[size(alig_seqs,2) size(alig_seqs,3)]);
                elseif size(alig_seqs,1) == 1
                    CM = reshape(alig_seqs,[size(alig_seqs,2) size(alig_seqs,3)]);
                else   % should not happen
                    error('Centroid is empty');
                end   
                if strcmp(version(end),'0')
                    %% Assign the centroids as the CM
                    Z{j} = CM;
                elseif strcmp(version(end),'1')
                    %% Assign the centroids as the closest sequence to the CM
                    costs = zeros(1,length(idx));
                    for i = 1:length(idx)
                        if strcmp(distance,'dtwCost')
                            W=dtwc(X{idx(i)},CM,1);
                            costs(i) = W(end,end);
                        elseif strcmp(distance,'euclidean')     % only possible if sequence lenghts are the same
                            costs(i) = sum(diag(pdist2(X{idx(i)},CM)))/length(X{idx(i)});                            
                        end
                    end
                    [~,closestCM] = min(costs);
                    Z{j} = X{idx(closestCM)};
                end
            end
            %% Convergence criterion: Centroids doesn't change
            if strcmp(distance,'dtwCost')
                WZ = dtwc(Z{j},Z2{j},1);
                costsZ(j) = WZ(end,end);
            elseif strcmp(distance,'euclidean')     % only possible if sequence lenghts are the same
                costsZ(j) = sum(diag(pdist2(Z{j},Z2{j})))/length(Z{j});                
            end
%             V(j,it+1) = costsZ(j);
            if costsZ(j) == 0 || isequal(Z{j},Z2{j}) || costsZ(j) == costs2Z(j)
                convergs(j) = true;
            else
                convergs(j) = false;
            end
            costs2Z(j) = costsZ(j);
%             text = strcat(text,sprintf(' %.2f',costsZ(j)));        
        end
%         fprintf('%s\n',text);
    end
    
    %% Expectation (assign new clusters)
    dc = inf*ones(k,length(X));
    for j = 1:k
        if isempty(Z{j})
            continue;
        end
        for i = 1:length(X)
            if isequal(X{i},Z{j})
                dc(j,i) = 0;
            else
                if strcmp(distance,'dtwCost')
                    W=dtwc(X{i},Z{j},1);
                    dc(j,i) = W(end,end);
                elseif strcmp(distance,'euclidean')     % only possible if sequence lenghts are the same
                    dc(j,i) = sum(diag(pdist2(X{i},Z{j})))/length(X{i});
                end
            end
        end
    end
    C = zeros(1,length(X));    
    for i = 1:length(X)
        [~,kf] = min(dc(:,i));
        C(i) = kf;
%         CG(i,it+1) = kf;
    end    
    
    it = it + 1;  
end

