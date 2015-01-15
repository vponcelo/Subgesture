function [X,Y] = readSkelsSeqs(seqs,nframesSeg,k)
% Read skeleton sequences and obtain the skeletons and labels for the 
% nubmer of given sequences
%
% output
%     X: cell with training samples at {1} and validation samples at {2}
%     Y: cell with training labels at {1} and validation labels at {2}
% input
%     seqs: list of skeleton sequences to consider
%     nframesSeg: fixed number of frames to segment the sequences '0' means no fixed segmentation ; '1' means frame segmentation
%     k: number of initial data clusters to visualize skeletons (if k > 0)

global VISUALIZE;
global DATATYPE;

    %% Read and save Skeleton data
%     if ~exist('data/',DATATYPE,'/skeletons.mat','file');
        paths{1} = strcat('C:/victor/thesis/datasets/',DATATYPE,'/samplesTrain/');
%         paths{1} = strcat('data/',DATATYPE,'/samplesTrain/');
        filenames{1} = dir(strcat(paths{1},'*.zip'));
        paths{2} = strcat('C:/victor/thesis/datasets/',DATATYPE,'/samplesVal/');
%         paths{2} = strcat('data/',DATATYPE,'/samplesVal/');
        filenames{2} = dir(strcat(paths{2},'*.zip'));
        paths{3} = strcat('C:/victor/thesis/datasets/',DATATYPE,'/samplesTest/');
%         paths{3} = strcat('data/',DATATYPE,'/samplesTest/');
        filenames{3} = dir(strcat(paths{3},'*.zip'));
        if ~sum(seqs)
            seqs = zeros(1,length(filenames));
            for i = 1:length(filenames)
                seqs(i) = length(filenames{i});
            end
        end
        videoData = cell(1,length(seqs));  % cell(1,2) without test data
        skeletons = cell(1,length(seqs));  % cell(1,2) without test data    
        
        display('Generating skeleton sequences...');
        for j = 1:length(videoData)
            videoData{j} = cell(1,seqs(j));
            skeletons{j} = cell(1,seqs(j));  
            for i = 1:length(videoData{j})
                samplePath = filenames{j}(i).name;
                fprintf('Generating %s ...\n',samplePath(1:end-4));
                if strcmp(DATATYPE,'chalearn2014')
%                     videoData{j}{i} = getChalearnData(strcat(paths{j},samplePath));
                    videoData{j}{i} = getChalearnDataNatalia(strcat(paths{j},samplePath));
                else
                    videoData{j}{i} = getSkeleton(strcat(paths{j},samplePath));
                end
                skeletons{j}{i}.skeleton = videoData{j}{i}.Skeleton;
                skeletons{j}{i}.angles = videoData{j}{i}.BufH;
                %skeletons{j}{i}.angles = videoData{j}{i}.WorldRotation;
                skeletons{j}{i}.labels = videoData{j}{i}.Labels;                
            end
        end        

    [features,skeletons] = getSkeletonFeatures( skeletons );
    X = cell(1,length(seqs));      % cell(1,2) without test data
    Y = cell(1,length(seqs));      % cell(1,2) without test data
    for i = 1:length(features)
%         X{i} = features{i};
        [X{i},Y{i}] = setLabels(features{i},skeletons{i},nframesSeg);        
    end
    
    %% Perform clustering over several training sequences for visualization
    if VISUALIZE;
        for i = 1:length(X)
            if k > 0,
                [~,~] = clusterSkels( X{1}, k, 800, skeletons{1}, nframesSeg );
            end
        end
    end
    
