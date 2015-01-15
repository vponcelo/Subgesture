function [in,fi]=detectSeqC(costMatrix)
% Get detected gesture with minimum costly path
% NOTE: The cost matrix should have dimension (n+1)x(m+1)
%    First row and column starts with 0 at position (1,1) and 
%    inifities at the remaining positions for aligning two sequences, and
%    0's on the first row when detecting gestures.

% input: 
%   costMatrix: dtw cost matrix
% output:
%   in: initial frame from the optimal path
%   fi: end frame from the optimal path

% 
if costMatrix(1,1) == 0 && unique(costMatrix(2:end,1) == inf)
    costMatrix(:,1) = [];    
else
    error('First column of the cost matrix must be 0 at the position (1,1) and infinity at the remaining positions (:,1)');
end

se = detectSeq_c(costMatrix);
in = se(1);
fi = se(2);
