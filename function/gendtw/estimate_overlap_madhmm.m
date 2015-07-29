function [O1,O2,R] = estimate_overlap_madhmm(GT,PREDSEQ,perc)

%% O1 - Without considering the iddle gesture as class
%% O2 - Considering the iddle gesture as class
%% perc - threshold on overlap to be considered as recognized

%% caso binario, 
O1=NaN;
O2=NaN;
% cGTtestkFr =GTtestkFr;
% cdetSeqLog=detSeqLog;
    
% cGTtestkFr(ofin)=1;
% cdetSeqLog(ofin2)=1;

ovs = sum(GT & PREDSEQ)./sum(GT | PREDSEQ);


ofin=find(GT==1);
ofin2=find(PREDSEQ==1);

R.rec=(length(intersect(ofin,ofin2)))./length(ofin);
R.prec=(length(intersect(ofin,ofin2)))./length(ofin2);


% clear ofin*
% ofin=find(SGT(:,3)==clx(i));
% ofin2=find(SPRED(:,3)==clx(i));
%     
% %% recall estimation
% O2c=zeros(1,length(ofin));
% for j=1:length(ofin),
%     cgest=cGTtestkFr(SGT(ofin(j),1):SGT(ofin(j),2));
%     pgest=cdetSeqLog(SGT(ofin(j),1):SGT(ofin(j),2));
%         inte=sum((cgest & pgest));
%         if (inte)./length(cgest) > perc,
%             O2c(j)=1;
%         else 
%             O2c(j)=0;
%         end                        
%     end
%     O2Pc=zeros(1,length(ofin2));
%     %% precision estimation    
%     for j=1:length(ofin2),
%         cgest=cGTtestkFr(SPRED(ofin2(j),1):SPRED(ofin2(j),2));
%         pgest=cdetSeqLog(SPRED(ofin2(j),1):SPRED(ofin2(j),2));
%         inte=sum((cgest & pgest));
%         if (inte)./length(cgest) > perc,
%             O2Pc(j)=1;
%         else 
%             O2Pc(j)=0;
%         end                        
%     end
%     Rec(i)=sum(O2c)./length(O2c);
%     Prec(i)=sum(O2Pc)./length(O2Pc);
%     
%     
% 
% % ovs(i) = sum(cGTtestkFr & cdetSeqLog)./sum(cGTtestkFr | cdetSeqLog);   
% 
% 
% fpth=110; %% number of frames to consider a gesture as detected (false pos.)
% clx=unique(GT); clx(clx==0)=[];
% O1=zeros(1,length(clx)-1);
% O2=zeros(1,length(clx));
% DETT=zeros(1,length(clx));
% 
% % GTtestkFr=zeros(size(GT));
% % detSeqLog=zeros(size(GT));
% 
% 
% for i=1:length(clx),
%     
%     
% %     overlapok
%     
%     ofin=find(GT==clx(i));
%     ofin2=find(PREDSEQ==clx(i));
%     
% %     cGTtestkFr =GTtestkFr;
% %     cdetSeqLog=detSeqLog;
% %     
% %     cGTtestkFr(ofin)=1;
% %     cdetSeqLog(ofin2)=1;
% %     
% %     ovs(i) = sum(cGTtestkFr & cdetSeqLog)./sum(cGTtestkFr | cdetSeqLog);   
% 
%     
%     nofin4=setdiff(1:length(GT),ofin);
%     
%     inter2=intersect(nofin4,ofin2);
%     DETT(i)=round(length(inter2)./fpth);
%         
%     inte=intersect(ofin,ofin2);
% 
%     if length(inte)./length(ofin) > perc,
%         O2(i)=1;
%     end
% end
% 
% gesture=find(clx~=-1);
% O1=O2(gesture);
% 
% R.prec=sum(O1)./(sum(O1)+sum(DETT(gesture)));
% R.prec2=sum(O2)./(sum(O2)+sum(DETT));
% 
% R.rec=sum(O1)./(length(clx)-1);
% R.rec2=sum(O2)./(length(clx));
% 
% % overlap=mean(ovs);
end