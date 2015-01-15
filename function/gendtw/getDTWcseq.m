function [s,in,fi] = getDTWcseq(W)
%%% This is a faster version to avoid multiple calls to detectSeq
%%% function.
% getDTWcseq returns the updated minimum costs given a DTW matrix
    % input:
    %   W: Warping cost sequence
    % output    
    %   s: array of updated minimum costs


%% Get costs Subsequences from DTW path
if numel(W) > 1
    if unique(W(1,2:end)) == 0 || unique(W(1,2:end)) == Inf
        W(1,:) = [];
    end
    if unique(W(2:end,1)) == Inf
        W(:,1) = [];
    end    
    s = inf*ones(1,size(W,2));
    
    % Forward step
    for fi = 1:length(s)
        current = size(W,1);
        in = fi;
        first = false;
        while current > 1 && ~first
            if in == 1
                if W(current-1,in) <= W(current,in)
                    first = true;
                end
            else
                if W(current-1,in-1) <= W(current-1,in) 
                    if W(current-1,in-1) <= W(current,in-1)
                        in = in-1;
                        current = current-1;
                    else
                        in = in-1;
                    end
                else 
                    if W(current-1,in) <= W(current,in-1)
                        current = current-1;
                    else
                        in = in-1;
                    end
                end
            end
        end
        
        % Intra-forward (review) step
%         fprintf('%d ',in);
        for i = in:fi
            if i==385036
                disp('');
            end
            if W(end,fi) < s(i)
                s(i) = W(end,fi);
            end
        end        
    end
%     if length(s) <= 4000,
%         disp(strcat(strrep(['Cost sequence: (' sprintf(' %.2f ', s ) ')'], ',)', ')')));       
%     end
%     disp('-------');
else
    error('getDTWcseq:costMatrix','matrix must be of dimension > 1');
end



