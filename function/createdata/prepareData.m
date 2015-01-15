function [X,Y,Xtest,Ytest] = prepareData(nrsamples,nseqs,nframesSeg,k)
% Generate the learning and test data from a set of options
%
% output
%     X: whole training data 
%     Y: training data labels
%     Xtest: Test data (empty means it has posterior assignment)
%     Ytest: test data labels (empty means it has posterior assignment)
% input
%     nrsamples: number of random samples for generating the random data
%     nseqs: number of skeleton sequences (or batches) to consider
%     nframesSeg: number of fixed frames per sequence segmentation. '0' means no fixed segmentation ; '1' means no subgestures (1 frame per subgesture).
%     k: number of data clusters (to visualize skeletons (if k > 0))

global DATATYPE;
global COORDS;
global PERCENTDATA;
global NORMTYPE;
global JOINTS;
global NAT;

switch DATATYPE
    case 'random'
        if ~exist('data/randomData.mat','file'),
            display('Generating data...');
            % Create random data
            dim=209;
            X=cell(1,nrsamples); 
            Y=cell(1,nrsamples); 
            % Create random states (3) and modify data to acomplish the following:
            %    State 1: Only 0-99
            %    State 2: Only 100-199
            %    State 3: Only 1000-1099
            for i=1:nrsamples,
                % Set the length for this sequence
                len=50+round(10*rand(1));

                % Generate the descriptors for this sequence
                seqData=round(rand(len,dim)*100);

                % Set the start and end frames
                v1=3+round(rand(1)*len*0.2);
                v2=3+round(rand(1)*len*0.2);

                % Modify descriptors
                %seqData(1:v1)=mod(seqData(1:v1)*2,10);
                %seqData(v1+1:end-v2-1)=mod(seqData(v1+1:end-v2-1),3)*3;
                %seqData(end-v2:end)=mod(seqData(end-v2:end)*2+1,10);
                seqData(1:v1)=mod(seqData(1:v1),100);
                seqData(v1+1:end-v2-1)=mod(seqData(v1+1:end-v2-1),100)+10;
                seqData(end-v2:end)=mod(seqData(end-v2:end),100)+1000;

                % Modify states
                seqStates=ones(1,len)*2;
                seqStates(1:v1)=1;
                seqStates(end-v2:end)=3;
                X{i}=seqData;
                Y{i}=seqStates;            
            end
            Xtest = X(nrsamples-nrsamples*0.2+1:nrsamples); X(nrsamples-nrsamples*0.2+1:nrsamples)=[];
            Ytest = Y(nrsamples-nrsamples*0.2+1:nrsamples); Y(nrsamples-nrsamples*0.2+1:nrsamples)=[];
            save('data/randomData.mat','dim','nsamples','X','Y','Xtest','Ytest');
            display('Done!');
        else
            display('Loading data...');
            load('data/randomData.mat');
        end
        display('Done!');
    
    case 'chalearn2013'
        display('Loading Chalearn 2013 data...');
        if ~exist(strcat('data/chalearn2013/chalearnData_',COORDS,'_',num2str(length(JOINTS)),'joints.mat'),'file'),
            [X,Y] = readSkelsSeqs(nseqs,nframesSeg,k);            
%             minValue = min(min(toMat(X)));
%             [maxValues,~] = max(toMat(X));
%             maxValue = max(maxValues);
            for i = 1:length(X)
%                 X{i} = (X{i}-minValue) / (maxValue-minValue);   % Data normalization
                Y{i} = setGT(X{i},Y{i},k,nframesSeg);
            end
            Xtest = X{3};   X(3) = [];
            Ytest = Y{3};   Y(3) = [];
            save(strcat('data/chalearn2013/chalearnData_',COORDS,'_',num2str(length(JOINTS)),'joints.mat'),'X','Y','Xtest','Ytest');
        else            
            load(strcat('data/chalearn2013/chalearnData_',COORDS,'_',num2str(length(JOINTS)),'joints.mat'));
%             [idxRowsDel,X] = cleanData(X);
%             if strcmp(COORDS,'world')
%                 load(strcat('data/chalearn2013/chalearnData_world_original','.mat'));
%                 for i = 1:length(idxRowsDel)
%                     X{i}(idxRowsDel{i},:) = 0;
%                 end
%                 save(strcat('data/chalearn2013/chalearnData_',COORDS,'_',num2str(length(JOINTS)),'joints.mat'),'X','Y','Xtest','Ytest');
%             end            
            if strcmp(NORMTYPE,'neck')
                X = norm2neck(X);
            elseif strcmp(NORMTYPE,'xyzangles')
                
            end
        end
        if PERCENTDATA < 100
            for i = 1:length(Y)
                ns = round(length(Y{i}.L)*PERCENTDATA/100);
                Y{i}.cnames = Y{i}.cnames(1:ns);
                Y{i}.seg = Y{i}.seg(1:ns+1);
                Y{i}.L = Y{i}.L(1:ns);
            end
        end
        display('Done!');
   case 'chalearn2014'
        switch NAT
            case 0
                s = 'joints.mat';
            case 1
                s = 'joints_natH.mat';
            case 2
                s = 'joints_nat.mat';
            case 3
                s = 'joints_natdepe3.mat';
        end
        display('Loading Chalearn 2014 data...');
        if ~exist(strcat('data/chalearn2014/chalearnData_',COORDS,'_',num2str(length(JOINTS)),s),'file'),
            [X,Y] = readSkelsSeqs(nseqs,nframesSeg,k);            
            for i = 1:length(X)
                Y{i} = setGT(X{i},Y{i},k,nframesSeg);
            end
            Xtest = X{3};
            Ytest = Y{3};
            save(strcat('data/chalearn2014/chalearnData_',COORDS,'_',num2str(length(JOINTS)),s),'X','Y','Xtest','Ytest');
        else            
            load(strcat('data/chalearn2014/chalearnData_',COORDS,'_',num2str(length(JOINTS)),s));
%              [idxRowsDel,X] = cleanData(X);
%             if strcmp(COORDS,'world')
%                 load(strcat('data/chalearn2014/chalearnData_',COORDS,'_',num2str(length(JOINTS)),'joints.mat'));
%                 for i = 1:length(idxRowsDel)
%                     X{i}(idxRowsDel{i},:) = 0;
%                 end
%                 save(strcat('data/chalearn2014/chalearnData_',COORDS,'_',num2str(length(JOINTS)),'joints.mat'),'X','Y','Xtest','Ytest');
%             end
            if strcmp(NORMTYPE,'neck')
                X = norm2neck(X);
                
                %% Hugo 
%                 miv=min(X{1});
%                 mav=max(X{1});
%                 mev=mav-miv;
% %                 X{1}(:,length(JOINTS)+1:end)=(X{1}(:,length(JOINTS)+1:end)-repmat(miv(length(JOINTS)+1:end),size(X{1},1),1))./(repmat(mev(length(JOINTS)+1:end),size(X{1},1),1));
%                 X{1}=(X{1}-repmat(miv,size(X{1},1),1))./(repmat(mev,size(X{1},1),1));
%                 miv=min(X{2});
%                 mav=max(X{2});
%                 mev=mav-miv;                
% %                 X{2}(:,length(JOINTS)+1:end)=(X{2}(:,length(JOINTS)+1:end)-repmat(miv(length(JOINTS)+1:end),size(X{2},1),1))./(repmat(mev(length(JOINTS)+1:end),size(X{2},1),1));
%                 X{2}=(X{2}-repmat(miv,size(X{2},1),1))./(repmat(mev,size(X{2},1),1));
            elseif strcmp(NORMTYPE,'xyzangles')
                
            end
        end
        if PERCENTDATA < 100
            for i = 1:length(Y)
                ns = round(length(Y{i}.L)*PERCENTDATA/100);
                Y{i}.cnames = Y{i}.cnames(1:ns);
                Y{i}.seg = Y{i}.seg(1:ns+1);
                Y{i}.L = Y{i}.L(1:ns);
            end
        end        
        display('Done!');
end