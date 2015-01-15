function addLabels(dataPath,labelsPath,dstPath)
    %
    % Add labels to unlabeled data. This script opens each ZIP file in the
    % data path and modify it SampleXXXXX_data.mat file in order to include
    % the labels, and create a new ZIP with updated data.
    %
    %    dataPath: Path of data files (it contains the files
    %              SampleXXXXX.zip)
    %
    %    labelsPath: Path of label files (it contains the files
    %              SampleXXXXX.txt)
    %    
    %    [dstPath:] This parameter is optional. If it is provided, updated
    %              ZIP files are stored in the given path, and source ZIP
    %              files remain unchanged. Otherwise, source ZIP files are
    %              modified.
    %
    
    
    % Check the parameters
    if nargin<2,
        error('Incorrect number of parameters');
    end
    
    % Prepare the paths    
    if exist('dataPath','var'),
        dataPath= regexprep(dataPath, '\\', '/');
        if dataPath(end)~='/',
            dataPath = strcat(dataPath,'/');
        end
    end
    if exist('labelsPath','var'),
        labelsPath= regexprep(labelsPath, '\\', '/');
        if labelsPath(end)~='/',
            labelsPath = strcat(labelsPath,'/');
        end
    end
    
    % Validate input paths
    if ~exist(dataPath,'dir'),
        error('Given dataPath must be a valid directory');
    end
    if ~exist(labelsPath,'dir'),
        error('Given labelsPath must be a valid directory');
    end
    
    % Prepare output path
    if ~exist('dstPath','var'),
        dstPath=dataPath;
    else
        dstPath= regexprep(dstPath, '\\', '/');
        if dstPath(end)~='/',
            dstPath = strcat(dstPath,'/');
        end
        if ~exist(dstPath,'dir'),
            mkdir(dstPath);
        end
    end
          
    % Read the files in the data directory
    files=dir(fullfile(dataPath,'Sample*.zip'));
    
    % Extract the labels from all the files
    for f=files',
        % Get the sample id
        [~,sampleID]=fileparts(char(f.name));
        
        % Read the labels file
        labFile=fullfile(labelsPath,[sampleID '.txt']);
        labels=getLabels(labFile);
        if isempty(labels),
            warning(['Cannot find labels for file ' sampleID]);
            continue;
        end       
        
        % Unzip the file
        dstZipFolder=sprintf('%s%s',dstPath,sampleID);
        recStat=recycle('off');
        if exist(dstZipFolder,'dir'),
            rmdir(dstZipFolder,'s');
        end
        unzip(fullfile(dataPath,f.name),dstZipFolder);
        
        % Read the data file (Contains Video structure)
        dataFile=fullfile(dstZipFolder,[sampleID '_data.mat']);
        
        if exist('Video','var'),
            clear Video;
        end
        load(dataFile);
        if ~exist('Video','var'),
            warning(['Error loading data file for sample' sampleID]);
            continue;
        else    
            % Store the labels
            Video.Labels=labels;
            save(dataFile,'Video');
        
            % Create the new zip
            zip(sprintf('%s.zip',dstZipFolder),'*.*',dstZipFolder);        
        end
        
        % Remove unziped folder
        rmdir(dstZipFolder,'s');
        recycle(recStat);
        
    end    
        
end


% Read the labels in the given file
function labels=getLabels(file),
    % Open the file
    fid=fopen(file,'r');
    if fid<0,
        labels=[];
        return;
    end
    
    % Read the file
    labels=[];
    while ~feof(fid),
        % Read next line
        line=fgets(fid);
        if line==-1,
            continue;
        end
                
        % Parse the line
        data=textscan(line,'%s','delimiter',';');
        data=data{:};
        
        % Store the information
        newLab.Name=char(data{1});
        newLab.Begin=str2num(data{2});
        newLab.End=str2num(data{3});
        labels=[labels newLab];
        
    end        
    
    % Close the file
    fclose(fid);

end


