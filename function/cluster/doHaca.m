function doHaca(X,Y)

% Input - 
%   X: Frame Data
%   Y: data labels

global FRAMECOMPRESS;

% The value of tag could be set to any integer between 1 and 14
% which correspods to the trial number of subject 86.
% You can derive similar results as shown in http://humansensing.cs.cmu.edu/projects/aca_more_results.html
tag = 14; 

%% source
wsSrc = mocSegSrc(tag);
[para, paraH] = stFld(wsSrc, 'para', 'paraH');

%% feature
if iscell(X),
    len = 40;
    x = []; y = [];
    if length(X) == length(Y)
        for i = 1:length(X)
            for j = 1:length(X{i})/len
                x = [x; X{i}(j*len-len+1:j*len)];
                y = [y; Y{i}(j)];
            end
        end
        X = x; Y = y; X = X'; Y = Y'; clear x y;
    end    
    s = 1;
    for c = Y.L,
        Indices = find(Y == c);
        s = [s Indices(end)];
    end
    clear Y c;
    segT = newSeg('s', s, 'G', L2G(Y.L, length(unique(Y.L))));
    cnames = {'gesture1','gesture2','gesture3','gesture4'};
else
    X = X';
    segT = newSeg('s', Y.seg, 'G', L2G(Y.L, length(unique(Y.L))));
    cnames = Y.cnames;
end
para.nMi = 6;
para.nMa = 10;
paraH(1).nMi = 2;
paraH(1).nMa = 4;
paraH(2).nMi = 2;
paraH(2).nMa = 4;
if FRAMECOMPRESS,
    wsData = mocSegData(wsSrc,X,segT,cnames);
    [X, segT] = stFld(wsData, 'X', 'segT');
end
K = conKnl(conDist(X, X), 'nei', .02);
para.nIni = 1;
para.k = length(unique(Y.L));
paraH(2).k = length(unique(Y.L));
paraH(1).k = paraH(2).k * 3;

%% init
seg0s = segIni(K, para);

%% aca
tic;
[segAca, segAcas] = segAlg('aca', [], K, para, seg0s, segT);
toc;

%% haca
tic;
seg0s = segIni(K, paraH(1));
segHaca = segAlg('haca', [], K, paraH, seg0s, segT);
toc;

%% gmm
tic;
segGmm = 0;
segGmm = segAlg('gmm', [], K, para, seg0s, segT);
toc;

%% save results
tic;
% save('/results/chalearn2012/visualClusters/haca/haca.mat','para','paraH','seg0s','X','K','segT','segAca','segHaca','segGmm');
toc;

clear K X;
display('Done!');

%% visualize results
showSegmentation(K,X,segT,segAca,segHaca,segGmm);

