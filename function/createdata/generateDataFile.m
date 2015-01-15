function stats=generateDataFile(outFile,trainPath,validationPath,testPath),
    close all;
    outFile='/Users/xbaro/Documents/Chalearn/val_tmp';
    trainPath='/Users/xbaro/Documents/Chalearn/Validation/';
    %validationPath='../Data/Validation/';
    %testPath='../Data/Test/';
    
    % Return% Usage(Training, PublicTest, PrivateTest
    
    % Prepare the paths    
    if exist('trainPath','var'),
        trainPath= regexprep(trainPath, '\\', '/');
        if trainPath(end)~='/',
            trainPath = strcat(trainPath,'/');
        end
    end
    
    
    if exist('validationPath','var'),
        validationPath= regexprep(validationPath, '\\', '/');
        if validationPath(end)~='/',
            validationPath = strcat(validationPath,'/');
        end
    end
    
    if exist('testPath','var'),
        testPath= regexprep(testPath, '\\', '/');
        if testPath(end)~='/',
            testPath = strcat(testPath,'/');
        end
    end
    
    % Open output file
    fout=fopen(outFile,'w');
       
    % Read training samples
    if exist('trainPath','var'),
        
        % Read the files in the directory
        files=rdir(trainPath,'*.zip');
        
        stats.trainning.numSamples=0;
        stats.trainning.numGestures=zeros(1,20);
        stats.numFiles=length(files);
        
    
        % Add sample data in the output file
        h = waitbar(0,'Exporting trainning examples');
        count=0;
        
        %procFiles={};
        %if exist('tmp.mat','file'),
        %    load('tmp.mat');
        %end
        
        for f=files',    
            %if ismember(f.name,procFiles),
            %    continue;
            %end
            tic;
            % Get the sample ID and labels for this sample
            [sampleID,labels,sampleStats]=getSampleLabels(fullfile(trainPath,f.name),1,0);
            
            % Accumulate the stats
            stats.trainning.numGestures=stats.trainning.numGestures+sampleStats;
            
            % Write data to output file
            writeSample(fout,sampleID,labels,'Training'); 
            
            % Add file to list
            %procFiles=[procFiles f.name];
            
            %save('tmp.mat','procFiles');
            count=count+1;
            waitbar(count/length(files),h,sprintf('Exporting trainning examples [%d/%d]',count,length(files)));
            toc
        end    
        close(h);
    end
    
    % Read validation samples
    if exist('validationPath','var'),
        
        % Read the files in the directory
        files=dir(strcat(validationPath,'*.zip'));
    
        % Add sample data in the output file
        for f=files',
            % Get the sample ID and labels for this sample
            [sampleID,labels]=getSampleLabels(fullfile(validationPath,f.name));
            
            % Write data to output file
            writeSample(fout,sampleID,labels,'PublicTest');            
        end    
    end
    
    % Read Test samples
    if exist('testPath','var'),
        
        % Read the files in the directory
        files=dir(strcat(testPath,'*.zip'));
    
        % Add sample data in the output file
        for f=files',
            % Get the sample ID and labels for this sample
            [sampleID,labels]=getSampleLabels(fullfile(testPath,f.name));
            
            % Write data to output file
            writeSample(fout,sampleID,labels,'PrivateTest');            
        end    
    end
    
    % Close the output file
    fclose(fout);
    
end

%%
% Get the labels seqüence for a given sample.
%
%     usage:
%         [id,labels]=getSampleLabels('../Data/Sample00001.zip');
%
function [sampleID,labels,stats]=getSampleLabels(samplePath,exportSkeleton,exportFrames)

    % If the video file is provided, open the video
    if ~exist('exportSkeleton','var'),
        exportSkeleton=0;
    end
    if ~exist('exportFrames','var'),
        exportFrames=0;
    end
   
    % Split the path and the filename
    [path,sampleID,ext]=fileparts(samplePath);
    
    %display(sampleID)
    % Check that the given file is a zip file
    if strcmp(ext,'.zip')~=1,
        warning('Expected a ZIP file as input');
    end
       
    % Check if there is a folder for this sample (unziped) or we need to unzip the file.
    sampleDir=fullfile(path,sampleID);
    src='zip';
    if exist(sampleDir,'dir'),
        src='dir';
    end
    
    % If it is necessary, unzip the file
    if strcmp(src,'zip'),
        try
           files=unzip(samplePath,sampleDir);
           if isempty(files),
               warning(['Empty extracted file: ' samplePath]);
           end
        catch
            warning(['Cannot extract file: ' samplePath]);
        end
            
    end
    
    % Read the data file
    data=load(fullfile(sampleDir,[sampleID '_data.mat']));
    
    % Read video
    if exportFrames==1,
        videoFile = VideoReader(fullfile(sampleDir,[sampleID '_color.mp4']));
    end
    
    
    % Use loaded data to generate the outout labels
    
    % Define an initial label for each frame (0-> no gesture)
    labels=zeros(1,data.Video.NumFrames);
    
    % Initialize the stats
    stats=zeros(1,20);
        
    % For each gesture in the sample generate corresponding labels
    for lab=data.Video.Labels,
        error=0;
        % Check the labels 
        if lab.Begin<=0 || lab.End<=0,
            if error==0,
                error=1;
                warning(['Null index errors in sample: ' samplePath ' (Gesture -> ' lab.Name ')']);
            end
            continue;
        end
        % Check that the index are correct
        if lab.Begin<1 || lab.End<1 || lab.End<=lab.Begin || lab.End>data.Video.NumFrames,
            warning(['Index error in sample: ' samplePath ' (Gesture -> ' lab.Name ')']);
            tmp=lab.Begin;
            lab.Begin=lab.End;
            lab.End=tmp;
        end
        % Check that all afected frames are 0
        if sum(labels(lab.Begin:lab.End))~=0,
            warning(['There are overlapped gestures in sample: ' samplePath]);
        end
        
        % Get the gesture identifier        
        id=getGestureID(lab.Name);
        
        % Add this value
        stats(id)=stats(id)+1;
                
        % Check that the gesture name is correct
        if id<=0,
            warning(['Unrecognized gesture: ' samplePath]);
        end
        
        % Set the labels
        labels(lab.Begin:lab.End)=id;  
        
        % Export the frames
        if exportFrames,
            exportPathRGB=sprintf('framesSample/%02d_%s',id,lab.Name);
            if ~exist(exportPathRGB,'dir'),
                mkdir(exportPathRGB);
            end
            exportGestureFrames(sprintf('%s/%s',exportPathRGB,sampleID),videoFile,lab);           
        end
        
        if exportSkeleton,
            exportPathSkel=sprintf('framesSample/%02d_%s',id,lab.Name);
            if ~exist(exportPathSkel,'dir'),
                mkdir(exportPathSkel);
            end
            exportSkeletonFrames(sprintf('%s/%s',exportPathSkel,sampleID),data.Video.Frames,lab);
        end
    end
    
    % If we unziped the file, remove the folder
    if strcmp(src,'zip'),
        recycleStat=recycle('off');
        try
            rmdir(sampleDir,'s');        
        catch err
            warning(['Cannot remove foler: ' sampleDir ' error: ' err.message]);
        end
       
        recycle(recycleStat);
    end
    
end
function exportSkeletonFrames(imageFile,skeleton,lab)
    numFrames=lab.End-lab.Begin;
    f=round(lab.Begin:numFrames/5:lab.End);
    s=0.5;
    bgImage=ones(480*s,640*s)*255;
    showSkeleton(skeleton(f(1)).Skeleton.PixelPosition*s,bgImage,1);
    f1=getframe;
    showSkeleton(skeleton(f(2)).Skeleton.PixelPosition*s,bgImage,1);
    f2=getframe;
    showSkeleton(skeleton(f(3)).Skeleton.PixelPosition*s,bgImage,1);
    f3=getframe;
    showSkeleton(skeleton(f(4)).Skeleton.PixelPosition*s,bgImage,1);
    f4=getframe;
    showSkeleton(skeleton(f(5)).Skeleton.PixelPosition*s,bgImage,1);
    f5=getframe;
    showSkeleton(skeleton(f(6)).Skeleton.PixelPosition*s,bgImage,1);
    f6=getframe;
  
    fout=[f1.cdata f2.cdata f3.cdata;f4.cdata f5.cdata f6.cdata];
    imwrite(fout,sprintf('%s_%d_%d.jpg',imageFile,lab.Begin,lab.End));
end
function exportGestureFrames(imageFile,videoFile,lab)
    numFrames=lab.End-lab.Begin;
    f=round(lab.Begin:numFrames/5:lab.End);
    s=0.5;
    f1=imresize(read(videoFile, f(1)),s);
    f2=imresize(read(videoFile, f(2)),s);
    f3=imresize(read(videoFile, f(3)),s);
    f4=imresize(read(videoFile, f(4)),s);
    f5=imresize(read(videoFile, f(5)),s);
    f6=imresize(read(videoFile, f(6)),s);
    fout=[f1 f2 f3;f4 f5 f6];
    imwrite(fout,sprintf('%s_%d_%d.jpg',imageFile,lab.Begin,lab.End));
end


%% 
% Write the sample line to the output file with the correct format
function writeSample(file,sampleID,labels,split)
    % Add the sample ID
    fprintf(file,'%s',sampleID);
    
    % Add the labels
    for l=labels,
        fprintf(file,',%d',l);
    end
    
    % Add the split
    
    % Add final line character
    fprintf(file,'\n');
end


%%
% Get all the files in a folder and subfolders of this folder
%
function files=rdir(path,pattern)

    % If no pattern is provided, assume all files
    if ~exist('pattern','var'),
        pattern='*.*';
    end
    
    % Prepare the path
    path= regexprep(path, '\\', '/');
    if path(end)~='/',
        path = strcat(path,'/');
    end
    
    % Get all the directories in this folder
    allFiles=dir(path);
    dirs=[];
    for d=allFiles',
        if d.isdir,
            dirs=[dirs;d];
        end
    end
    
    % Get all the files in root folder
    files=[];
    fdir=dir(sprintf('%s%s',path,pattern));
    for file=fdir',
        files=[files;file];
    end
    
    % Add the files in the subfolders
    for d=allFiles',
        if strcmp(d.name,'.')==1 || strcmp(d.name,'..')==1,
            continue;
        end
        dfiles=dir(sprintf('%s%s/%s',path,d.name,pattern));
        for file=dfiles',
            file.name=sprintf('%s/%s',d.name,file.name);
            files=[files;file];
        end
    end
    

end
%%
% Get the gesture identifier from their string description
%
%   usage:
%       id=getGestureID('vieniqui');
function gestureID=getGestureID(gestureName)
    
    % Define the list of available gestures in order (1....N) 
    gestureNameList={'vattene','vieniqui','perfetto','furbo','cheduepalle','chevuoi','daccordo','seipazzo','combinato','freganiente', ...
                     'ok','cosatifarei','basta','prendere','noncenepiu','fame','tantotempo','buonissimo','messidaccordo','sonostufo'};
            
    % Get the position in the array of gestures
    gestureID=find(strcmp(gestureNameList,gestureName));
    
    % Check that this gesture exists
    if isempty(gestureID),
        warning(['Cannot identify gesture: ' gestrureName]);
    end
    
end