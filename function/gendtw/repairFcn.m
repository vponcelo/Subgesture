function RP = repairFcn(P, params, maxSegs, nvars)
% function that receives a population and repairs the invalid genes with
% brute force, given parameters.

RP = P;

for i = 1:size(RP,1)
    idxi = logical([0 mod(2:size(RP,2),2) == 0]);                       % start seq indices
    idxf = logical([0 mod(2:size(RP,2),2) == 1]);                       % length seq indices
    
    %% Fix bad initial-length segments
    badIni = RP(i,:) < 1 | (RP(i,:) > maxSegs-params.nmax & RP(i,:) < inf); 
    badIni = badIni & idxi;                                                 % Detect out of bounds initial segments
    RP(i,badIni) = randi([params.k0 maxSegs-params.nmax+1],1,sum(badIni));  % force new start of bad initial segments
    lensExceed = find(RP(i,:) > maxSegs & RP(i,:) < inf & idxf);            % Indices that exceed the maximum X length
    if ~isempty(lensExceed)
        RP(i,lensExceed(1)-1) = maxSegs-params.nmax;                        % force new start of length exceeded segments
        RP(i,lensExceed(1)) = params.nmax;                                  % force new length of length exceeded segments
        if length(lensExceed) > 1
            RP(i,lensExceed(2:end)-1) = inf; RP(i,lensExceed(2:end)) = inf;
        end
    end
    lensExceed = find((RP(i,:) < params.nmin | RP(i,:) > params.nmax & RP(i,:) < inf) & idxf);    % Indices that exceed the maximum segment length
    if ~isempty(lensExceed)
        RP(i,lensExceed) = randi([params.nmin params.nmax],1,length(lensExceed));                    % force new length of length exceeded segments        
    end
    
    %% Fix k if it's out of bounds
    k = RP(i,1); kf = find(RP(i,:)==inf,1)-1;
    if ~(k >= params.k0 && k <= nvars-1), k = randi([params.k0 nvars-1]); end   % force new k 
    while ~(k >= params.k0 && k <= kf),
        if (k-params.k0)/(kf-params.k0) > 1                                     % k trends to kf
            if k-5 > params.k0                                                  
                k = k - 5;
            elseif k-5 < params.k0
                k = k - 1;
            end
        else
            if kf + 10 < nvars                                                  % k trends to k0
                RP(i,kf+1:2:kf+10) = randi([1 maxSegs-params.nmax+1],1,5);      % force new start point
                RP(i,kf+2:2:kf+10) = randi([params.nmin params.nmax],1,5);      % force new length
                kf = kf + 10;
            else
                RP(i,kf+1) = randi([1 length(X)-params.nmax+1]);                % force new start point
                RP(i,kf+2) = randi([params.nmin params.nmax]);                  % force new length
                kf = kf + 2;
            end
        end
    end
    RP(i,1) = k;
end

end