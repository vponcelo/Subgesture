function W = dtwc(x,t,align,w,D,KM,KT)
% returns the DTW matrix

% input
%    x: vector you are testing against
%    t: vector you are testing
%    align: flag indicating if DTW is used to align two secuences
%       (align=1) or to detect start-end gestures (align=0)
%    w: Window size
%    D: Dissimilarity matrix to use instead of another distance metric
%    KM: Cost matrix of representing the model sequence into subgestures
%    KM: Cost matrix of representing the test sequence into subgestures
    
% output:
    % D: warping cost matrix
 
%% Dynamic Time Warping
if nargin == 3 
    W=dtw_c(double(x),double(t),logical(align));
elseif nargin == 4
    W=dtw_c(double(x),double(t),logical(align),w);
elseif nargin == 7
    W=dtw_c(double(x),double(t),logical(align),w,double(D),double(KM),double(KT));
else
    error('Number of arguments is incorrect;');
end
% it = 0;
% while sum(sum(W(2:end,2:end) < 0)) > 0 && it < 5
%     W=dtw_c(x,t,align);
%     it = it + 1;
%     if it >= 5
%         W = dtw2(x,t,logical(align));
%     end
% end
W(:,1) = inf;
if align
    W(1,:) = inf;
end
W(1,1) = 0;

% W2=[x; W'];               % For 1-D vectors only
% visualDTW=[[0 t]' W2];    % For 1-D vectors only
