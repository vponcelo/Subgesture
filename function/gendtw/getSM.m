function M = getSM(params,Xtrain_l,model,mnseg,mk)
% function that obtain the Subgesture Models

display('Computing Subgesture Models for each gesture...');
if strcmp(params.mType,'modelSM1')
    % Median/Max Models to new Subgestures Models
    if nargin < 5
        error('getSM:nArgs','Number of arguments must be at lesst 5');
    end
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
            warning('getSM:kIncorrect','number of clusters cannot be greater than number of segments. Decreasing k (not consistent).');
        end
        [~,~,mErrsV,~,timeV,~,SM] = runKMeansDTW(params.version,mk,'dtwCost',mk,[],[],[],[],[],Msegs,[]);
        [~,kV] = min(mErrsV);
        emptyCells = cellfun(@isempty,SM{kV}{timeV});
        SM{kV}{timeV}(emptyCells) = [];
        M{ng} = SM{kV}{timeV}{end};
    end    
elseif strcmp(params.mType,'allSM1')
    % Gesture samples to new Median Subgesture Models
    if nargin < 5
        error('getSM:nArgs','Number of arguments must be at lesst 5');
    end    
    SM = cell(1,length(params.M));
    for ng = 1:length(params.M)
        SM{ng} = cell(1,length(Xtrain_l{ng}));
        for ns = 1:length(Xtrain_l{ng})
            XngSegs = cell(1,mnseg);
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
            warning('getSM:kOutOfBound','number of clusters cannot be greater than number of segments. Decreasing k (not consistent).');
        end
        [~,~,mErrsV,~,timeV,~,Z] = runKMeansDTW(params.version,mk,'dtwCost',mk,[],[],[],[],[],XngSegs,[]);
        [~,kV] = min(mErrsV);
        emptyCells = cellfun(@isempty,Z{kV}{timeV});
        Z{kV}{timeV}(emptyCells) = [];
        SM{ng}{ns} = Z{kV}{timeV}{end};
    end
    M = getModels(SM,length(SM),params.mType,false,params.usemax_l);
elseif strcmp(params.mType,'modelSM2')
    % Median Models to the corresponding Median Subgestures Models
    M = cell(1,length(params.M));
    for i = 1:length(params.M)
        alig_seqs = zeros(length(model.SM),size(params.M{i},1),size(params.M{i},2));
        for j = 1:length(model.SM)
            if params.resize,
                alig_seqs(j,:,:)=imresize(model.SM{j},[size(params.M{i},1),size(model.SM{j},2)]);
            else
                W = dtwc(model.SM{j},params.M{i},1);
                [~,~,alig_seqs(j,:,:)]=aligngesture(model.SM{j},W);
            end
        end
        M{i} = reshape(mean(alig_seqs),[size(alig_seqs,2) size(alig_seqs,3)]); 
    end
elseif strcmp(params.mType,'allSM2')
    % Gesture samples to the corresponding Median Subgesture Models
    M = cell(1,length(params.M));
    for i = 1:length(params.M)
        M{i} = cell(1,length(Xtrain_l{i}));
        display(sprintf('Aligning samples of gesture %d to Subgestures',i));
        for sg = 1:length(Xtrain_l{i})
            alig_seqs = zeros(length(model.SM),size(Xtrain_l{i}{sg},1),size(Xtrain_l{i}{sg},2));
            for j = 1:length(model.SM)
                if params.resize,
                    alig_seqs(j,:,:)=imresize(model.SM{j},[size(Xtrain_l{i}{sg},1),size(model.SM{j},2)]);
                else
                    W = dtwc(model.SM{j},Xtrain_l{i}{sg},1);
                    [~,~,alig_seqs(j,:,:)]=aligngesture(model.SM{j},W);
                end
            end
            M{i}{sg} = reshape(mean(alig_seqs),[size(alig_seqs,2) size(alig_seqs,3)]);
        end
    end
    M = getModels(M,length(M),params);
else
    error('getSM:OptionError','Option chosen for the Subgesture Models is incorrect. Check the value of params.mType');
end