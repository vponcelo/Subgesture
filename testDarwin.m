function s = testDarwin(X,Xt,Y,Yt)

%% parameters
CVAL=1;
clop_model=kridge;

if iscell(X) && iscell(Xt)
    display('Processing training ...')
    k = 1;
    for i = 1:length(X)
        for j = 1:length(X{i})
            W = genRepresentation(X{i}{j}',CVAL);
            if mod(i,10)==0, display(sprintf('Current gesture %d/%d of class %d/%d\n',j,length(X{i}),i,length(X))); end
            Xtrain(k,:)=W';
            Ytrain(k,1)=i;
            k = k + 1;
        end
    end
    display('Processing test ...')
    k = 1;
    for i = 1:length(Xt)
        for j = 1:length(Xt{i})
            W = genRepresentation(Xt{i}{j}',CVAL);
            if mod(i,10)==0, display(sprintf('Current gesture %d/%d of class %d/%d\n',j,length(Xt{i}),i,length(Xt))); end
            Xtest(k,:)=W';
            Ytest(k,1)=i;
            k = k + 1;
        end
    end
else
    ntr=length(Y.seg);
    nte=length(Yt.seg);
    display('Processing training ...')
    for i=2:ntr,
        if i == ntr, fi = Y.seg(i); else fi = Y.seg(i)-1; end
        X1=X(Y.seg(i-1):fi,:);
        W = genRepresentation(X1,CVAL);
        if mod(i,10)==0, display(sprintf('Current gesture segment %d of %d\n',i,ntr)); end
        Xtrain(i-1,:)=W';
        Ytrain(i-1,1)=Y.L(i-1);
    end
    display('Processing test ...')
    for i=2:nte,
        if i == nte, fi = Yt.seg(i); else fi = Yt.seg(i)-1; end
        X1=Xt(Yt.seg(i-1):fi,:);    
        W = genRepresentation(X1,CVAL);
        if mod(i,10)==0, display(sprintf('Current gesture segment %d of %d\n',i,nte)); end
        Xtest(i-1,:)=W';
        Ytest(i-1,1)=Yt.L(i-1);
    end
end
%% classfication part
classex=unique(Ytrain);
nclassex=length(classex);
Yones=-ones(size(Ytrain,1),nclassex);
YTones=-ones(size(Ytest,1),nclassex);
for i=1:nclassex,
    Yones(Ytrain==classex(i),i)=1;
    YTones(Ytest==classex(i),i)=1;
end
[a,b]=train(one_vs_rest(clop_model),data(Xtrain,Yones));
[c]=test(b,data(Xtest,YTones))
pred=zeros(size(Ytest));
for i=1:nclassex,
    pred(c.X(:,i)==1)=classex(i);
end

s = (length(find((pred-Ytest)==0))./length(Ytest));

% % 51694300 / corp. santander, 5485
