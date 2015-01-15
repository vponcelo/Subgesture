% Copyright (C) 2013 Víctor Ponce-López <vponcel@uoc.edu>,
% Scene UNderstanding and Artificial Intelligence Laboratory,
% Internet Interdisciplinary Institute (IN3)
% Universitat Oberta de Catalunya, 08018 Barcelona, Spain.

function [in,fi,alignedgesture]=aligngesture(gesture,costMatrix)
% Get aligned gesture with the minimum cost 
% NOTE: The cost matrix should have dimension (n+1)x(m+1)
%    First row and column starts with 0 at position (1,1) and
%    inifities at the remaining positions when aligning two sequences, or
%    0's on the first row when detecting gestures.

% input: 
%   gesture: reference sequence for the alignment 
%       (use an empty list to detect start-end of gesture)
%       
%   costMatrix: dtw cost matrix
% output:
%   in: initial frame from the optimal path
%   fi: end frame from the optimal path
%   alignedgesture: aligned sequence of minimum cost 
%       (returns an empty list for start-end detection)

% 
%% initialization
if costMatrix(1,1) == 0 && unique(costMatrix(2:end,1) == inf)
    costMatrix(:,1) = [];
    if unique(costMatrix(1,:)) == inf
        costMatrix(1,:) = [];
    end
else
    error('First column of the cost matrix must be 0 at the position (1,1) and infinity at the remaining positions (:,1)');
end

if costMatrix(1,:) >= 0
    minpos = cell(1,size(costMatrix,1));
    k = size(costMatrix,1);
    minpos{k} = size(costMatrix,2);
    fi = minpos{k}; % set end frame
else
    error('Values of the cost matrix cannot be negative');
end
i = size(costMatrix,1); j = size(costMatrix,2);
if j == 1
    for pos = 1:length(minpos)
        minpos{pos} = 1;
    end
end

%% begin backward search of optimal the path 
setPos = false;
pos = zeros(1,j);
while isempty(minpos{1})    
    if i > 1 && j > 1
        M = costMatrix(i-1:i,j-1:j);
        M = [M(1,1) M(1,2) M(2,1)];    
        [~,minp] = min(M);
        switch minp            
            case 1  % diagonal backward
                j = j - 1;
                i = i - 1;
                setPos = true;                
            case 2  % row backward
                i = i - 1;
                setPos = true;
            case 3  % column backward
                j = j - 1;                
        end
        
        % assign positions
        if i > 1
            pos(j) = j;
            if i == 2 && isempty(gesture)   % ~any(costMatrix(1,:))
                i = 1;
                setPos = true;
                minpos(i) = [];
                k = k - 1;
            end
        end
    end
    
    % first row reached
    if i == 1 && j > 1 && ~isempty(gesture) % any(costMatrix(1,:))
        j = j - 1;
        pos(j) = j;
        setPos = false;
        % first column reached
        if j == 1
            setPos = true;
        end
    end
    
    % set positions
    if setPos       
        if j == 1 && (i > 1 || ~any(pos))
            if i >= k
                error('Cost Matrix is incorrect');
            end
            pos = j;
            i = i - 1;
        else
           setPos = false;
        end
        k = k - 1;
        pos(pos==0) = [];
        minpos{k} = pos;
        pos = zeros(1,j);
    end
end
in = minpos{1}(end);    % set start frame

%% Use empty gesture to detect start-end
if ~isempty(gesture)    
    alignedgesture = zeros(length(minpos),size(gesture,2));
    %% Begin alignment by computing the means on each aligned frame
    for i = 1:length(minpos)
        if length(minpos{i}) > 1
            alignedgesture(i,:) = mean(gesture(minpos{i},:));
        else
            alignedgesture(i,:) = gesture(minpos{i},:);
        end
    end
else
    alignedgesture = [];
end