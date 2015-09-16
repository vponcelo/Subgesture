function [D,U] = computeSGFromU(Us,X,Y)
% input:
    % Us: Subgesture sequences of U
    % X: sequence or set of sequences for each class
    % Y: data labels
% output:
    % D: Subgesture representations of X in Us
    % U: Weighted Matrix of the updated ranking machine

if iscell(X)
    D = cell(1,length(X));
    U = cell(1,length(X));
    for m = 1:length(X)
        D{m} = cell(1,length(X{m}));
        U{m} = cell(1,length(X{m}));
%         m
        for j = 1:length(X{m})
            U{m}{j} = inf*ones(size(Us,1),size(X{m}{j},1));
            sws = 10:10:size(X{m}{j},1);
            for i = 1:length(sws)
                for st = 1:sws(i):size(X{m}{j},1)-1
                    fi = min(st+sws(i)-1,size(X{m}{j},1));
                    if fi > 1
                        seq = X{m}{j}(st:fi,:);
                        W = genRepresentation(seq,1);
                        for s = 1:size(Us,1)
                            d = pdist2(W',Us(s,:));
                            for c = st:fi
                                if d < U{m}{j}(s,c), U{m}{j}(s,c) = d; end
                            end
                        end
                    end
                end
            end
            [~,D{m}{j}] = min(U{m}{j}); 
        end
    end
else
    if size(X,1) > 10000 && ~exist('Y','var')
        Xs = X; t = 500;        % t fixed parts. Other possibility is to discretize each gesture from its label
        if mod(length(Xs),t) > 0,
            ns = round(length(Xs)/t) + 1;
        else
            ns = length(Xs)/t;
        end
    elseif exist('Y','var') 
        Xs = X; ns = length(Y.L); offset = 1;
    else
        ns = 1;
    end
    U = cell(1,ns); 
    for n = 1:ns
        if ns > 1,
            if ~exist('Y','var')
                X = Xs(1+(n-1)*t:min(n*t,length(Xs)),:); 
            else
                if n == length(Y.L), offset = 0; end
                X = Xs(Y.seg(n):Y.seg(n+1)-offset,:);
            end
        end
        U{n} = inf*ones(size(Us,1),size(X,1));
        sws = 10:10:length(X);
        for i = 1:length(sws)
            for st = 1:sws(i):size(X,1)-1
                fi = min(st+sws(i)-1,size(X,1));
                if fi > 1
                    seq = X(st:fi,:);
                    W = genRepresentation(seq,1);
                    for s = 1:size(Us,1)
                        d = pdist2(W',Us(s,:));
                        for c = st:fi
                            if d < U{n}(s,c), U{n}(s,c) = d; end
                        end
                    end
                end
            end
        end
    end
    emptyCells = cellfun(@isempty,U); U(emptyCells) = [];
    U = cell2mat(U);
    [~,D] = min(U);
end