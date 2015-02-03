function M = getMSM(params,Xtrain_l,mnseg,mk)
% function that obtain the Median Subgesture Model

display('Computing Median Subgesture Models for each gesture...');
if strcmp(params.mType,'directMSM1')
    % Median Models to Median Subgesture Models
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
        if mk > length(Msegs)
            error('getMSM:ManySegs','number of clusters cannot be greater than number of segments');
        end
        [~,~,mErrsV,~,timeV,~,MSM] = runKMeansDTW(params.version,mk,'dtwCost',mk,[],[],[],[],[],Msegs,[]);
        [~,kV] = min(mErrsV);
        emptyCells = cellfun(@isempty,MSM{kV}{timeV});
        MSM{kV}{timeV}(emptyCells) = [];
        M{ng} = MSM{kV}{timeV}{end};
    end    
elseif strcmp(params.mType,'directMSM2')
    % Gesture samples to Median Subgesture Models
    MSM = cell(1,length(params.M));
    for ng = 1:length(params.M)
        MSM{ng} = cell(1,length(Xtrain_l{ng}));
        for ns = 1:length(Xtrain_l{ng})
            if strcmp(params.msmType,'fix')
                XngSegs = cell(1,mnseg);
                idx = 1; fseg = round(size(Xtrain_l{ng}{ns},1)/mnseg); % fixed segmentation
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
            if mk > length(XngSegs)
                error('getMSM:ManySegs','number of clusters cannot be greater than number of segments');
            end
            [~,~,mErrsV,~,timeV,~,Z] = runKMeansDTW(params.version,mk,'dtwCost',mk,[],[],[],[],[],XngSegs,[]);
            [~,kV] = min(mErrsV);
            emptyCells = cellfun(@isempty,Z{kV}{timeV});
            Z{kV}{timeV}(emptyCells) = [];
            MSM{ng}{ns} = Z{kV}{timeV}{end};
        end                
    end
    M = getMedianModels(MSM,length(MSM),params.mType,false);
end
display('Done!');