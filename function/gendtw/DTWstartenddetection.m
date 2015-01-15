%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GMM-DTW Implementation
%Usage -> [start_point,end_point,cost,M]=DTW(sec,ptr,cost,dist,align)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%       sec: the secuence in which the pattern is searched (feautures are columwise)
%       ptr: the searched pattern (feautures are columwise), for the case
%       of 'probability' ptr is a cell in which every position is a
%       gmdistribution object.
%       cost: the cost threshold
%       dist: 'euclidean','probability','histogram' if distance is
%       'probability' then sec is a two column matrix with the mean (column 1)
%       and variance (column 2) and as many rows as frames
%       align: flag indicating if DTW is used to align two secuences
%       (align=1) or to detect start-end gestures (align=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output:
%       start_point: starting point of the pattern
%       end_point: end point of the pattern
%       cost: cost of the particular pattern
%       THE ABOVE THREE OUPUTS ARE REPEATED FOR AS MANY PATTERNS FOUND
%       WITH LESS COST THAN THE INPUT COST THRESHOLD
%       M: cost matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [start_point,end_point,cost,M]=DTWstartenddetection(sec,ptr,cost,distance,align)

%% Initialize M matrix accordingly
if align
    M=zeros(size(ptr,1),size(sec,1));
    M(1,:)=inf;
    M(:,1)=inf;
    M(1,1)=0;

   
else
    M=zeros(size(ptr,1),size(sec,1));
end

%% main DTW loop
for i=2:size(M,2)
    for j=2:size(M,1)
        
        if strcmp('euclidean',distance)
                dist_dtw=pdist2(sec(i,:),ptr(j,:),'euclidean');
        
        elseif strcmp('probability',distance)
            
                model_prob=mean(mahal(ptr{j},sec(i,:))-log(ptr{j}.PComponents)*2);
                dist_dtw=model_prob;

        elseif strcmp('histogram',distance)
            dist_dtw=hist_dist(sec(i,:),ptr(j,:));
            
        end
%         if i>1 && j>1
            M(j,i)=(dist_dtw+min([M(j-1,i),M(j,i-1),M(j-1,i-1)]));
%         else
%             M(j,i)=dist_dtw;
%         end
    end
end


%% Obtain warping paths in M 
start_point=backward_pth(M,find(M(end,:)<cost));
end_point=find(M(end,:)<cost);

%% Filter warping paths
[start_point,end_point]=filter_paths(start_point,end_point);
cost=M(end,end_point);


end


function start_point=backward_pth(M,points)
start_point=zeros(size(points,1));

for p=1:size(points,2)
    
    current=[size(M,1),points(p)];
    while current(1)>1 && current(2)>1
        
        [next,pos]=min([M(current(1)-1,current(2)),M(current(1),current(2)-1),M(current(1)-1,current(2)-1)]);
        if pos==1
            current=[current(1)-1,current(2)];
        elseif pos==2
            current=[current(1),current(2)-1];
        else
            current=[current(1)-1,current(2)-1];
        end
    end
    start_point(p)=current(2);
    
end

end

function dist=hist_dist(histA,histB)

dist=sum(min(histA,histB));
end

function [start_p,end_p]=filter_paths(start_point,end_point)

for i=2:length(start_point)
    if start_point(i)<=end_point(i-1)
        start_point(i)=-1;
        end_point(i-1)=-1;
    end
end
start_p=start_point(find(start_point>0));
end_p=end_point(find(end_point>0));
end

