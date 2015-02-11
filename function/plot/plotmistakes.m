function plotmistakes(predictions, Yval, saveflg)
% 
%   Plots the errors of a model. A segment of the whole sequence is shown. 
%
%   Inputs: 
%       * predictions:  1D vector of predictions
%       * Yval:         structure with ground truth labels.
%
%

% size of the segment to sho
seglen=1000;
% get predictions in a Nx3 matrix
P = transform_predictions(predictions);
% get ground truth labels
GTtestk=Yval.Lfr;
% discard sequences of length smaller than 10
gelen=P(:,2)-P(:,1);
ofin=find(gelen>10);
P=P(ofin,:);
% zeros variable for predctions
predica=zeros(1,length(GTtestk));
% randonmly select the segment to be shown
r=randperm(round(length(GTtestk)./2));
segment=r(1):r(1)+seglen;

% start figure
figure;gcf; hold on;
set(gca,'FontSize',14);

% plot results for each class
for hi=1:20,
    ofin=find(P(:,3)==hi);
    cpred=predica;
    gtc=predica;
    gtc(find(GTtestk==hi))=1;
    for hj=1:length(ofin),
        cpred(P(ofin(hj),1):P(ofin(hj),2))=1;
    end
    ovlpc(hi) = sum(gtc & cpred)/sum(gtc | cpred);
    subplot(4,5,hi);plot([gtc(segment);]','LineWidth',1.5); 
    gcf;hold on;
    plot([cpred(segment);]','.-r','LineWidth',1.5); ylim([-1 2]);
    ylim([-1 2]);xlim([0.5 seglen]);
    err=(abs([gtc(segment)-cpred(segment);]));
    plot((find(err==1)),ones(1,length(find(err==1))),'gx');    
    title(['Ovl: ' num2str(ovlpc(hi))]);
end
legend({'GT','Pred','Errors'});
set(gcf,'Color','w');

if saveflg,
saveas(gcf,strrep(strrep(datestr(now),':','-'),' ',''),'jpg');
end
end
