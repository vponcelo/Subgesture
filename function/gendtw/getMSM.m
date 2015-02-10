function M = getMSM(params,Xtrain_l,model,mnseg,mk)
% function that obtain the Max/Median Subgesture Models

display('Computing Median Subgesture Models for each gesture...');
if strcmp(params.msmType,'fix') && strcmp(params.mType,'modelMSM1')
    % Median Models to new Median Subgestures Models
    M = cell(1,length(params.M));
    for ng = 1:length(params.M)
        Msegs = cell(1,mnseg);
        idx = 1; fseg = round(size(params.M{ng},1)/mnseg); % fixed segmentation
        for nsgs = 1:mnseg
            if nsgs < mnseg
                Msegs{nsgs} = params.M{ng}(idx:nsgs*fseg,:);
            else
                Msegs{nsgs} = params.M{ng}(idx:end,:);
            end
            idx = nsgs*fseg + 1;
        end
        emptyCells = cellfun(@isempty,Msegs);
        Msegs(emptyCells) = [];
        while mk > length(Msegs) && mk >= params.k0-1
            mk = mk - 1;
            warning('getMSM:kIncorrect','number of clusters cannot be greater than number of segments. Decreasing k (not consistent).');
        end
        try
            [~,~,mErrsV,~,timeV,~,MSM] = runKMeansDTW(params.version,mk,'dtwCost',mk,[],[],[],[],[],Msegs,[]);
        catch e
            error('getMSM:kIncorrect',e.message);
        end
        [~,kV] = min(mErrsV);
        emptyCells = cellfun(@isempty,MSM{kV}{timeV});
        MSM{kV}{timeV}(emptyCells) = [];
        M{ng} = MSM{kV}{timeV}{end};
    end    
elseif strcmp(params.msmType,'fix') && strcmp(params.mType,'allMSM1')
    % Gesture samples to new Median Subgesture Models
    MSM = cell(1,length(params.M));
    for ng = 1:length(params.M)
        MSM{ng} = cell(1,length(Xtrain_l{ng}));
        for ns = 1:length(Xtrain_l{ng})
            if strcmp(params.msmType,'fix')
                try
                    XngSegs = cell(1,mnseg);
                catch e
                    display(e.message)
                end
                idx = 1; fseg = round(size(Xtrain_l{ng}{ns},1)/mnseg)-1; % fixed segmentation
                for nsgs = 1:mnseg
                    if nsgs < mnseg
                        XngSegs{nsgs} = Xtrain_l{ng}{ns}(idx:nsgs*fseg,:);                        
                    else
                        XngSegs{nsgs} = Xtrain_l{ng}{ns}(idx:end,:);
                    end
                    idx = nsgs*fseg + 1;
                end
            end
            emptyCells = cellfun(@isempty,XngSegs);
            XngSegs(emptyCells) = [];
            while mk > length(XngSegs) && mk >= params.k0-1
                mk = mk - 1;
                warning('getMSM:kOutOfBound','number of clusters cannot be greater than number of segments. Decreasing k (not consistent).');
            end
            try
                [~,~,mErrsV,~,timeV,~,Z] = runKMeansDTW(params.version,mk,'dtwCost',mk,[],[],[],[],[],XngSegs,[]);
            catch e
                error('getMSM:kIncorrect',e.message)
            end            
            [~,kV] = min(mErrsV);
            emptyCells = cellfun(@isempty,Z{kV}{timeV});
            Z{kV}{timeV}(emptyCells) = [];
            MSM{ng}{ns} = Z{kV}{timeV}{end};
        end                
    end
    M = getMedianModels(MSM,length(MSM),params.mType,false,params.usemax_l);
elseif strcmp(params.msmType,'evoSegs') && strcmp(params.mType,'modelMSM2')
    % Median Models to the corresponding Median Subgestures Models
    M = cell(1,length(params.M));
    for i = 1:length(params.M)
        alig_seqs = zeros(length(model.SM),size(params.M{i},1),size(params.M{i},2));
        for j = 1:length(model.SM)
            W = dtwc(model.SM{j},params.M{i},1);
            [~,~,alig_seqs(j,:,:)]=aligngesture(model.SM{j},W);                
        end
        M{i} = reshape(mean(alig_seqs),[size(alig_seqs,2) size(alig_seqs,3)]); 
    end
elseif strcmp(params.msmType,'evoSegs') && strcmp(params.mType,'allMSM2')
    % Gesture samples to the corresponding Median Subgesture Models
    M = cell(1,length(params.M));
    for i = 1:length(params.M)
        M{i} = cell(1,length(Xtrain_l{i}));
        for sg = 1:length(Xtrain_l{i})
            alig_seqs = zeros(length(model.SM),size(Xtrain_l{i}{sg},1),size(Xtrain_l{i}{sg},2));
            for j = 1:length(model.SM)
                W = dtwc(model.SM{j},Xtrain_l{i}{sg},1);
                [~,~,alig_seqs(j,:,:)]=aligngesture(model.SM{j},W);
            end
            M{i}{sg} = reshape(mean(alig_seqs),[size(alig_seqs,2) size(alig_seqs,3)]);
            % tengo celda de celdas de ptrs
        end
        M = getMedianModels(M,length(M),params.mType,false,params.usemax_l);
    end        
    
else
    error('getMSM:OptionError','Otion selected to create the Subgesture Models is incorrect. Check params.mType and params.msmType variables');
end
display('Done!');