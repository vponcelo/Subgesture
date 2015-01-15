function createSeqVideo(S,outputPath,s)
% Create a video sequence to show the skeletons over selected frames
% input:
%   S: Skeleton sequence: Concatenation of 20 x coordinates plus 20 y
%       coordinates. If S is a cell, k sequences are created where k =
%       length(S) (each sequence should be a centroid sequence).
%   outputPath: Output directory path
%   s: scaling factor (window size) 0 < s <= 1. Default: s=1.


if ~exist('outputPath','var')
    error('Output path must be specified');
end
if ~exist('s','var')
    s = 1;
elseif s <= 0 || s > 1
    error('The scaling factor s must be greater than 0 and lower or equal than 1');
end

bgImage=ones(480*s,640*s)*255;

if iscell(S)    
    for i = 1:length(S)       
        videofile = sprintf('%s/%d.avi',outputPath,i);                    
        if ~exist(videofile,'file')
            aviobj = avifile(videofile,'compression','None');
            aviobj.Fps = 20;
            for j = 1:size(S{i},1)
                Zdraw = [S{i}(j,1:end/2); S{i}(j,end/2+1:end)]'*size(bgImage,2);
                if Zdraw > 0
                    h=showSkeleton(Zdraw*s,bgImage,1);
                    aviobj = addframe(aviobj, h);
                    close(h);
                end
            end    
            aviobj = close(aviobj);
        end
    end
else
    aviobj = avifile(sprintf('%s/seq2test',outputPath),'compression','None');
    aviobj.Fps = 20;
    for i = 1:size(S,1)
        Zdraw = [S(i,1:end/2); S(i,end/2+1:end)]'*size(bgImage,2);
        if Zdraw > 0
            h=showSkeleton(Zdraw*s,bgImage,1);
            aviobj = addframe(aviobj, h);
            close(h);
        end
    end
    aviobj = close(aviobj);
end
