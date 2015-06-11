
if isnan(X{1}), display(sprintf('%d',sum(sum(isnan(X{1}))))); error('NaN Data in X1'); end
if isnan(X{2}), display(sprintf('%d',sum(sum(isnan(X{2}))))); error('NaN Data in X2'); end
if isinf(X{1}), display(sprintf('%d',sum(sum(isinf(X{1}))))); error('Inf Data in X1'); end
if isinf(X{2}), display(sprintf('%d',sum(sum(isinf(X{2}))))); error('NaN Data in X2'); end

% Y{1}.L(Y{1}.L == 1000 | Y{1}.L == 1001) = 36;
% Y{2}.L(Y{2}.L == 1000 | Y{2}.L == 1001) = 36;

for i = 1:2
    if ~isequal(unique(Y{1}.L),unique(Y{2}.L)), error('Error in Labels'); end
    if ~isequal(unique(Y{1}.cnames),unique(Y{2}.cnames)), error('Error in label names'); end
    if ~isequal(Y{i}.seg(end),length(X{i})), error('Inconsistency in the last segment of Y{%i}',i); end
    for j = 2:length(Y{i}.seg)
        if Y{i}.seg(j-1) >= Y{i}.seg(j), 
            error('Segment inconsistency in Y{%i}',i); 
        end
    end
end

e = inf; s = inf;
while s+e >= length(Y{1}.seg)
    s = randperm(length(Y{1}.seg),1);
    e = s+round(0.3*length(Y{1}.seg));
end

Xtest = X{2}; Ytest = Y{2}; X{2} = []; Y{2} = [];

Y{2}.seg = Y{1}.seg(s:e);
Y{2}.seg = Y{2}.seg-Y{2}.seg(1)+1;
Y{2}.seg(end) = Y{2}.seg(end)-1;
Y{2}.L = Y{1}.L(s:e-1);
Y{2}.cnames = Y{1}.cnames(s:e-1);

X{2} = X{1}(Y{1}.seg(s):Y{1}.seg(e)-1,:);
if ~isequal(size(X{2},1),Y{2}.seg(end)), error('Segment error'); end
X{1}(Y{1}.seg(s):Y{1}.seg(e)-1,:) = [];

segs = Y{1}.seg(e:end)-Y{1}.seg(e)+Y{1}.seg(s);
Y{1}.seg(s:end) = []; Y{1}.seg = [Y{1}.seg segs];
Y{1}.L(s:e-1) = [];
Y{1}.cnames(s:e-1) = [];

if ~isequal(unique(Y{1}.L),unique(Y{2}.L),unique(Ytest.L)) || ...
        ~isequal(unique(Y{1}.cnames),unique(Y{2}.cnames),unique(Ytest.cnames))
    error('Error in Labels'); 
end
if ~isequal(Y{1}.seg(end),length(X{1})), error('Inconsistency in the last segment of 1'); end
if ~isequal(Y{2}.seg(end),length(X{2})), error('Inconsistency in the last segment of 2'); end
if ~isequal(Ytest.seg(end),length(Xtest)), error('Inconsistency in the last segment of test'); end

% save(strcat('data/MAD/MAD_5'),'X','Y','Xtest','Ytest');