function m_dtw = getModels(X,k,params)
% Get the refrence models from K-Means DTW
% 
% Input:
%   X: data sequences
%   k: number of models (clusters)
%   C: array with the cluster indices for which each sequence of X belongs
%   to.
%   params: set of parameters
% Ouptut:
%   m_dtw: cell with means of the k clusters

m_dtw = cell(1,k);
if strcmp(params.mType,'DCSR')
    W_ini = cell(1,k); 
end

display(sprintf('Computing gesture Models from training for the m=%d distinct gestures ...',k));

% Compute the warping costs W for all possible combinations without repetition
for i = 1:k
    if params.k > 0
        m_dtw{i} = cell(1,params.k);
        if isempty(m_dtw{i})
            m_dtw{i} = [];
        else
            [m_dtw{i},~] = kmeansDTW(X{i},params.k,params.version,params.dist,[]);
        end
    elseif strcmp(params.mType,'DCSR')
        % compute the warping costs W among all disctinct combinations
        % without repetition, then compute the warping costs to the
        % mean warped sequence.
        idxW = 1;
        W_ini{i} = cell(1,factorial(length(X{i}))/(factorial(2)*factorial(length(X{i})-2)));
        for j = 1:length(X{i})
            for k = 1:length(X{i})
                if j < k
                    D = dtwc(X{i}{j},ptr,1);
                    idxW = idxW + 1;
                end
            end
        end
        slengths = zeros(1,length(W_ini{i}));        
        for j = 1:length(W_ini{i})
            slengths(j) = size(W_ini{i}{j},1);
        end
        if ~params.usemax_l,
            if mod(length(slengths),2) == 0
                meanLength = mean(slengths);
                dists2mean = abs(slengths-meanLength);
                [~,idx] = min(dists2mean);
            else
                medianLength = median(slengths);
                idx = slengths == medianLength;
            end
        else
            [~,idx]=max(slengths);
        end
        m_dtw = W_ini;
        % At this point a normalization and a dimensionality reduction
        % technique is required for all W_ini having the same dimensionality
        ptr = W_ini{i}{idx}; 
        for j = 1:length(W{i})
            if j ~= idx
                [~,~,~,W]=DTWstartenddetection(W_ini{i}{j},ptr,0,'euclidean',1);                             
            else
                W = ptr;
            end                    
        end
    else
        m_dtw{i} = cell(1,length(X{i}));
        if isempty(m_dtw{i})
            m_dtw{i} = [];
        else
            % compute the warping costs W between the secuences and the
            % median sequence
            if iscell(X{i})
                len = length(X{i});
                slengths = zeros(1,len);
                for j = 1:length(X{i})
                    slengths(j) = size(X{i}{j},1);
                end
            else
                len = length(X);
                slengths = zeros(1,len);
                for j = 1:length(X)
                    slengths(j) = size(X{j},1);
                end
            end
            if ~params.usemax_l,
                if mod(length(slengths),2) == 0
                    meanLength = mean(slengths);
                    dists2mean = abs(slengths-meanLength);
                    [~,idx] = min(dists2mean);
                else
                    medianLength = median(slengths);
                    idx = find(slengths == medianLength);
                end         
                if length(unique(idx)) > 1
                    idx = round(length(slengths)/2);
                end
            else
                 [~,idx]=max(slengths);
%                 [~,idx]=min(slengths);
            end
            ptr = X{i}{idx}; 
            alig_seqs = zeros(length(X{i}),size(ptr,1),size(ptr,2));
            for j = 1:length(X{i})
                if j ~= idx                  
                    if params.resize,
                        alig_seqs(j,:,:)=imresize(X{i}{j},[size(ptr,1),size(X{i}{j},2)]);
                    else
                        W = dtwc(X{i}{j},ptr,1);
                        [~,~,alig_seqs(j,:,:)]=aligngesture(X{i}{j},W);
                    end                    
                else
                   alig_seqs(j,:,:)=ptr;
                end            
            end        
            
            % Compute the mean among the aligned warping matrices
            if ~params.gmm
                if size(alig_seqs,1) > 1
                    m_dtw{i} = reshape(mean(alig_seqs),[size(alig_seqs,2) size(alig_seqs,3)]);
                elseif size(alig_seqs,1) == 1
                    m_dtw{i} = reshape(alig_seqs,[size(alig_seqs,2) size(alig_seqs,3)]);
                end                
            else
                % use mean GMM distribution instead of the mean 
                GMMs = zeros(length(X{i}),size(ptr,1),size(ptr,2));
                for j = 1:length(W{i})
                    W_fit = reshape(alig_seqs(j,:,:),[size(alig_seqs,2) size(alig_seqs,3)]);
                    GMMs(j,:,:) = gmdistribution.fit(W_fit,3);
                end
                m_dtw{i} = reshape(mean(GMMs),[size(GMMs,2) size(GMMs,3)]);
            end
        end
    end
end