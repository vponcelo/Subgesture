function K = getUpdatedCosts(X,SM)
% Compute dtw matrices and obtain the updated warping costs between mean 
% subsequence models and the test sequence 
% Input:
%   X: sequence to test
%   SM: subsequence models
% Output:
%   MW: Normalized warping costs for each k

if ~iscell(SM) 
    error('Subsequence models must be a cell array');
end

nSM = length(SM);
K = inf*ones(nSM,size(X,1));
% tic;
for i = 1:nSM
    if ~isempty(SM{i})
        W = single(dtwc(X,SM{i},false));
%         tic;
%         for j = 1:length(W)-1
%             [in,~] = detectSeqC(W(:,1:j+1));
% %             [in2,~,~] = aligngesture([],W(:,1:j+1));
% %             [~,in3,~] = getDTWcseq(W(:,1:j+1));
% %             if ~isequal(in,in2,in3)
% %                 error('asdf');
% %             end
%             fprintf('%d ',in);
%             for k = in:j
%                 if W(end,j+1) < K(i,k)
%                     K(i,k) = W(end,j+1);            
%                 end
%             end
%         end
%         toc;
        if size(W,1) == 1
            display('W calculated so that SM is not empty');
            error('prog:input','dim(X)=%dx%d and dim(SM)=%dx%d\n', ...
                size(X,1),size(X,2),size(SM{i},2),size(SM{i},2));            
        end
%         tic;
%         [K(i,:),~,~] = getDTWcseq(W);     % This works faster        
%         toc;
%         tic;
        K(i,:) = getDTWcseq_c(W);
%         toc;
%         if ~isequal(K(i,:),a)
%             if sum(K(i,:) ~= a) > 1
%                 %realmax('single')==a(328032)
%                 error('getUpdatedCosts:diffUpdates','More than one differences detected between update costs functions');
%             end
%         end
    else
        K(i,:) = inf;        
    end
end
% toc;
% display('');