function dataClusters = getDataClusters(numClusters,numIterations,X)
% kmlsample -d 209 -k 60 -s 100 -max 1000 -df descData.out -out fileOut.txt
    X = toMat(X);
    cd('bin');
    dlmwrite('descData.out',X,' ');
    cmd=sprintf('kmlsample.exe -d %d -k %d -s %d -max %d -df descData.out -out fileOut.txt',size(X,2),numClusters,numIterations,size(X,1));
    system(cmd);        % File not found?
    dataClusters = dlmread('bin/fileOut.txt',';');
    
    delete('fileOut.txt');
    delete('descData.out');
    cd('..');