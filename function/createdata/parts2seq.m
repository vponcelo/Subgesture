dirpath = '/home/cvc/hupba/Data/VideoDarwin/Hollywood2/feats/trj/';
partfilepath = '/home/cvc/hupba/Datasets/Hollywood2/train_test_split.mat';

display('Listing directory with mbh-fv descriptors ...');
filenames = dir(dirpath);
X = cell(1,2); Y = cell(1,2);
X{1} = []; X{2} = []; Xtest = [];
Y{1}.seg = 1; Y{2}.seg = 1; Ytest.seg = 1;
Y{1}.L = []; Y{1}.L = []; Ytest.L = [];

display('Loading train/test labels ...');
load(partfilepath)
classid = labels2; clear labels2;
trn_indx = classid(cur_train_indx{1},:); clear cur_train_indx;
test_indx = classid(cur_test_indx{1},:); clear cur_test_indx;

% k=200;
for i=3:length(filenames)
    display(sprintf('Loading descriptor %s  ...',filenames(i).name));
    load(strcat(dirpath,filenames(i).name));
    %% uncomment / comment to perform different setting
    data = histv; clear histv;
%     display('Computing PCA ...');
%     [coeff,score,latent] = princomp(histv,'Algorithm','eig','NumComponents',500);
%     percent_explained = 100*latent/sum(latent);
%     display(sprintf('Percentage of data explanation with PCA: %.2f',sum(percent_explained(1:k))));
%     data = score(:,1:k); 
    %%
    
    if strcmp(filenames(i).name(12),'r')
        %train
        X{1} = [X{1}; data];
        if i < length(filenames)
            Y{1}.seg = [Y{1}.seg Y{1}.seg+size(data,1)];
            Y{1}.L = [Y{1}.L find(trn_indx(i,:) > 0)];
        else
            Y{1}.seg = [Y{1}.seg Y{1}.seg+size(data,1)-1];
            Y{1}.L = [Y{1}.L find(trn_indx(i,:) > 0)];
        end
        display(sprintf('Current training sequence length: %i frames\n',Y{1}.seg(end)));
    elseif strcmp(filenames(i).name(12),'e')
        %test
        Xtest = [Xtest; data];
        if i < length(filenames)
            Ytest.seg = [Ytest.seg Ytest.seg+size(data,1)];
            Ytest.L = [Ytest.L find(test_indx(i,:) > 0)];
        else
            Ytest.seg = [Ytest.seg Ytest.seg+size(data,1)-1];
            Ytest.L = [Ytest.L find(test_indx(i,:) > 0)];
        end
        display(sprintf('Current test sequence length: %i frames\n',Ytest.seg(end)));
    else
        error('part2seq:fname','incorrect filename');
    end
end

% Ideal setting. Compute PCA over all observations together (requires much more memory demand)
% len = size(X{1},1);
% Xtemp = [X{1};Xtest]; X{1} = []; Xtest = [];
% display('Computing PCA ...');
% [coeff,score,latent] = princomp(Xtemp,'Algorithm','eig','NumComponents',500);
% percent_explained = 100*latent/sum(latent);
% display(sprintf('Percentage of data explanation with PCA: %.2f\n',sum(percent_explained(1:k))));
% data = score(:,1:k);
% X{1} = data(1:l,:); Xtest = data(l+1:end,:);
display('Computing LDA ...');
options = [];
% options.Fisherface = 1;
[eigvector, ~] = LDA(Y{1}.L, options, X{1});
X{1} = X{1}*eigvector;

% save hollywood2_mbh-fv_pca.mat X Y Xtest Ytest -v7.3
save hollywood2_mbh-fv_lda.mat X Y Xtest Ytest -v7.3
