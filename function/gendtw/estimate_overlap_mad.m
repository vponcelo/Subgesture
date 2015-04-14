function [O1,O2,R] = estimate_overlap_mad(GT,PREDSEQ,perc)

%% O1 - Without considering the iddle gesture as class
%% O1 - Considering the iddle gesture as class
%% perc - threshold on overlap to be considered as recognized

fpth=110; %% number of frames to consider a gesture as detected (false pos.)
clx=unique(GT); clx(clx==0)=[];
O1=zeros(1,length(clx)-1);
O2=zeros(1,length(clx));
DETT=zeros(1,length(clx));

for i=1:length(clx),
    
    ofin=find(GT==clx(i));
    ofin2=find(PREDSEQ==clx(i));
    
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

end