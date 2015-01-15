function [videoData]=getSkeleton(samplePath)
%%
% Get all the information for a given sample.
%
%     videoData: Information structure. Contains the fields:
%        .NumFrames: Number of frames
%        .FrameRate: 20
%        .MaxDepth: Maximum depth value
%        .Labels: Array of gestures information.
%        .RGB: Cell array with all the RGB frames
%        .Depth: Cell array with all the Depth frames
%        .UserSegmentation: Cell array with the segmentation of the user
%                           for each frame
%        .Skeleton: Skeleton structure for each frame
%        .Audio: Structure with the audio samples and frequency.
%     labels: List of gesture identifiers for given sequence
%
%
%     usage:
%         [data,labels]=getSampleData('Sample00001.zip');
%
    global JOINTS;
    
    % Split the path and the filename
    [path,sampleID,ext]=fileparts(samplePath);
    
    % Check that the given file is a zip file
    if strcmp(ext,'.zip')~=1,
        warning('Expected a ZIP file as input');
    end
       
    % Check if there is a folder for this sample (unziped) or we need to unzip the file.
    sampleDir=fullfile(path,sampleID);
    src='zip';
%     if exist(sampleDir,'dir'),
%         src='dir';
%     end
    
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
    
    % Get the video information
    videoData=load(fullfile(sampleDir,[sampleID '_data.mat']));
    videoData=videoData.Video;
    videoData.SampleID=sampleID;
    % Get the Skeleton data
    for i = 1:videoData.NumFrames
        if length(JOINTS) == 20
            videoData.Skeleton(i)=videoData.Frames(i).Skeleton;
        else
            videoData.Skeleton(i).WorldPosition = videoData.Frames(i).Skeleton.WorldPosition(JOINTS,:);
            videoData.Skeleton(i).WorldRotation = videoData.Frames(i).Skeleton.WorldRotation(JOINTS,:);
            videoData.Skeleton(i).JointType = videoData.Frames(i).Skeleton.JointType(JOINTS);
            videoData.Skeleton(i).PixelPosition = videoData.Frames(i).Skeleton.PixelPosition(JOINTS,:);
            
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
    
    % showSkeleton(skeleton(#frame).Skeleton.PixelPosition*s,bgImage,1);
    % getframe;
    
end