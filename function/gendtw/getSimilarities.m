function S = getSimilarities(X)
% compute normalized similarity matrix from the distances between sequences
%   Input -
%       X: Cell with sequences to compare
%   Output -
%       S: Similary matrix
n = length(X);
A = zeros(n,n);
if n > 1
    for i = 1:n
        for j = 1:length(X)
            if i < j
                W = dtwc(X{i},X{j},1);
                A(i,j) = W(end,end);
            end 
        end
    end
    B = A./max(max(A));
    S = B + B';
else
    S = 1;
end


% imagesc(S)
% hold on;
% colormap(gray)
% hold off;
