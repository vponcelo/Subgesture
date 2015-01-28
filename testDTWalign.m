X = rand(80,100);       % test sequence
T = rand(43,100);       % pattern sequence

W = dtw(X,T,false);     % compute DTW cost matrix W, which is used afterwards for start-end detection
[in,fi,~] = aligngesture([],W);         % detect start-end (if not specified, fi is always the last position of W)

W = dtw(X,T,true);      % compute DTW cost matrix W, which is used afterwards for sequence alignment
[~,~,alignedgesture] = aligngesture(gesture,W); 