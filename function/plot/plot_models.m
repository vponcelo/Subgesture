function [h] = plot_models(M)
% draws the models, assuming one per gesture
%
% Input:  
%       M - a cell with each element being a prototype model
% Output:
%       h -  handler of the figure
%
%% plot the models
nmodels=length(M);
[A]=factor(nmodels);
d1=prod(A(1:end-1));
d2=(A(end));
cnt=1;
h=figure;
gcf;hold on;
for i=1:d1,
    for j=1:d2,
        subplot(d1,d2,cnt);imagesc(M{cnt});
        modsize(cnt)=size(M{cnt},1);
        cnt=cnt+1;
        
    end
end
%set(h,'Title','Models..');%title('Models...');

%% estimate the distance matrix (like a confusion matrix)
% D = dtwc(X{i}{j},ptr,1);
[~,ms]=max(modsize);
ptr=M{ms};
alig_seqs = zeros(nmodels,size(ptr,1),size(ptr,2));
for j = 1:nmodels
    if j ~= ms                   
    	W = dtwc(M{j},ptr,1);
        [~,~,alig_seqs(j,:,:)]=aligngesture(M{j},W);
    else
        alig_seqs(j,:,:)=ptr;
    end            
end

for i = 1:nmodels,
    uno=reshape(alig_seqs(i,:,:),[size(alig_seqs,2),size(alig_seqs,3)]);
    for j = 1:nmodels
        dos=reshape(alig_seqs(j,:,:),[size(alig_seqs,2),size(alig_seqs,3)]);
        Dis(i,j)=sum(sum(abs(uno-dos)));
    end
end
figure; imagesc(Dis);colorbar;
gcf; hold on; title('Estimated distance among models')
end