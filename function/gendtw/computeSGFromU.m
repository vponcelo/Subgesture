function [D,U] = computeSGFromU(Us,X)
% input:
    % Us: Subgesture sequences of U
    % X: sequence or set of sequences for each class
    % params: set of parameters
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
    U = inf*ones(size(Us,1),size(X,1));
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
                        if d < U(s,c), U(s,c) = d; end
                    end
                end
            end
        end
    end
    [~,D] = min(U);
end