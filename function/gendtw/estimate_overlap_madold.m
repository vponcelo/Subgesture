function [O1,O2,R] = estimate_overlap_madold(GT,PREDSEQ,perc)

%% O1 - Without considering the iddle gesture as class
%% O2 - Considering the iddle gesture as class
%% perc - threshold on overlap to be considered as recognized

fpth=110; %% number of frames to consider a gesture as detected (false pos.)
clx=unique(GT); clx(clx==0)=[];
O1=zeros(1,length(clx)-1);
O2=zeros(1,length(clx));
DETT=zeros(1,length(clx));

% GTtestkFr=zeros(size(GT));
% detSeqLog=zeros(size(GT));


for i=1:length(clx),
    
    
%     overlapok
    
    ofin=find(GT==clx(i));
    ofin2=find(PREDSEQ==clx(i));
    
%     cGTtestkFr =GTtestkFr;
%     cdetSeqLog=detSeqLog;
%     
%     cGTtestkFr(ofin)=1;
%     cdetSeqLog(ofin2)=1;
%     
%     ovs(i) = sum(cGTtestkFr & cdetSeqLog)./sum(cGTtestkFr | cdetSeqLog);   

    
    nofin4=setdiff(1:length(GT),ofin);
    
    inter2=intersect(nofin4,ofin2);
    DETT(i)=round(length(inter2)./fpth);
        
    inte=intersect(ofin,ofin2);

    if length(inte)./length(ofin) > perc,
        O2(i)=1;
    end
end

gesture=find(clx~=-1);
O1=O2(gesture);

R.prec=sum(O1)./(sum(O1)+sum(DETT(gesture)));
R.prec2=sum(O2)./(sum(O2)+sum(DETT));

R.rec=sum(O1)./(length(clx)-1);
R.rec2=sum(O2)./(length(clx));

% overlap=mean(ovs);
end