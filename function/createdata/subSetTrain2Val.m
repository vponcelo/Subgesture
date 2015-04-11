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
    error('Segment error'); 
end

% save(strcat('data/MSRACT3D/MSRACT3D_PCA'),'X','Y','Xtest','Ytest');