function S = findsegments(Seq)


ix=1;
cusym=Seq(ix);
S(1,[1:3])=0;
S(1,1)=ix;
S(1,2)=0;
S(1,3)=cusym;
nidx=1;
lseq=length(Seq);
while ix<lseq
    
    validx=ix+1:lseq;
    dSeg=find(Seq(validx)~=cusym);
    if ~isempty(dSeg),
        ix=validx(dSeg(1));
        cusym=Seq(ix);
        nidx=nidx+1;
        S(nidx,1)=ix;
        S(nidx,2)=0;
        S(nidx,3)=cusym;    
        S(nidx-1,2)=S(nidx,1)-1;
    
    else
        if ix<lseq,
            S(nidx,2)=lseq;
        end
        break;
    end
    
end

end