function [hmmTR,hmmE,hmmStates,pTrain,pVal] = learnModel(Dtrain,Ctrain,Xtrain,Xval,model)

% Output:
%   pTrain: estimated probabilities of each training sample
%   pVal: estimated probabilities of each test sample
%   hmmTR: estimated Transition Matrix
%   hmmE: estimated Emission Matrix,


% Input:
%   Dtrain: Discrete Training Data
%   Dtest: Discrete Test Data
%   Ctrain: Clustered Training Data
%   Ctest: Clustered Test Data
%   Xtrain: Training Data
%   Xval: Test Data

switch model
    case 1,
    % for s=1:3,
        hmmStates = 3;
        [hmmTR, hmmE]=learnHMM(hmmStates,Dtrain);
        pTrain = evaluateSequences(Ctrain,Xtrain,hmmTR,hmmE);
        pVal = evaluateSequences(Ctrain,Xval,hmmTR,hmmE);   
    % end
end

%% Create an S states HMM and estimate their probabilities
% if ~exist('hmmPars.mat','file'),
%     display('Estimating HMM parameters...');
%     [hmmTR, hmmE]=learnHMM(4,discreteData);
%     save('hmmPars.mat','hmmTR','hmmE');
%     display('Done!');
% else
%     load('hmmPars.mat');
% end
% 
% %% Apply HMM over training sequences
% if ~exist('resProbTrainSeqs.mat','file'),
%     display('Evauating HMM over training sequences...');
%     trainProbs=evaluateSequences(dataClusters,data,hmmTR, hmmE);
%     save('resProbTrainSeqs.mat','trainProbs');
%     display('Done!');
% else
%     load('resProbTrainSeqs.mat');
% end
% 
% %% Create test sequences
% if ~exist('testData.mat','file'),
%     display('Generating test data...');
%     % Create random data
%     numTest=100;
%     testData={};
%     testStates={};
%     % Create random states (3) and modify data to acomplish the following:
%     %    State 1: Only 0,2,4,6,8,10
%     %    State 2: Only 1,3,6,9
%     %    State 3: Only 1,3,5,7,9
%     for i=1:numTest,
%         % Set the length for this sequence
%         len=50+round(10*rand(1));
%         
%         % Create some initial well generated samples and others random
%         % samples
%         numCorrect=0.5;
%         if i<(numTest*numCorrect),
%             % Generate the descriptors for this sequence
%             seqData=round(rand(len,dim)*100);
%         
%             % Set the start and end frames
%             v1=3+round(rand(1)*len*0.2);
%             v2=3+round(rand(1)*len*0.2);
% 
%             % Modify descriptors
%             seqData(1:v1)=mod(seqData(1:v1),100);
%             seqData(v1+1:end-v2-1)=mod(seqData(v1+1:end-v2-1),100)+10;
%             seqData(end-v2:end)=mod(seqData(end-v2:end),100)+1000;
%             
%             % Modify states
%             seqStates=ones(1,len)*2;
%             seqStates(1:v1)=1;
%             seqStates(end-v2:end)=3;
%             
%             % Store this position
%             numTestCorrect=i;
%         else
%             % Generate the descriptors for this sequence
%             seqData=round(rand(len,dim)*1100);
%             seqStates=ones(1,len)*(-1);
%         end
%         
%         testData{i}=seqData;
%         testStates=seqStates;
%     end
%     save('testData.mat','numTest','testData','testStates','numTestCorrect');
%     display('Done!');
% else
%     load('testData.mat');
% end
% 
% %% Estimate the probabilities of test sequences
% if ~exist('resProbTestSeqs.mat','file'),
%     display('Evauating HMM over training sequences...');
%     testProbs=evaluateSequences(dataClusters,testData,hmmTR, hmmE);
%     save('resProbTestSeqs.mat','testProbs');
%     display('Done!');
% else
%     load('resProbTestSeqs.mat');
% end
% 
% %% Show results
% subplot(2,2,1);hist(trainProbs);title('Histogram of probabilities for training data');
% subplot(2,2,2);imagesc(hmmE);title('Labels distribution for each State');
% subplot(2,2,3);hist(testProbs(1:numTestCorrect));title('Histogram of probabilities for correct test data');
% subplot(2,2,4);hist(testProbs(numTestCorrect+1:end));title('Histogram of probabilities for random test data');
