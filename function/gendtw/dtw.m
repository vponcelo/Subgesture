% Computes the dynamic time warping matrix of two multi-dimensional signals

% Copyright (C) 2013 Víctor Ponce-López <vponcel@uoc.edu>,
% Scene UNderstanding and Artificial Intelligence Laboratory, IN3 (UOC),
% Dept. of Applied Mathematics and Analysis, University of Barcelona,
% Human Pose Recovery and Behavior Analysis, Computer Vision Center (UAB),
% Universitat Oberta de Catalunya, 08018 Barcelona, Spain.

% Original version for 1D signals from:
% Copyright (C) 2013 Quan Wang <wangq10@rpi.edu>,
% Signal Analysis and Machine Perception Laboratory,
% Department of Electrical, Computer, and Systems Engineering,
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

function W=dtw(s,t,a,w)
% Output:
    % W: DTW cost matrix
% Input:
    % s: test sequence
    % t: pattern sequence
    % a: flag either for alignment 'true' or sequence detection 'false'
    % w: window parameter
    %      if s(i) is matched with t(j) then |i-j|<=w
    % d: resulting distance

if nargin<4
    w=Inf;
    if nargin<3
        a=true;
    elseif ~islogical(a)
        error('dtw:param','a parameter must be a logical value');
    end
end

        % N-Dimensions
if size(s,2) > 1 && size(t,2) > 1
    ns=size(s,1);
    nt=size(t,1);
else    % 1-Dimension
    ns=length(s);
    nt=length(t);
end
w=max(w, abs(ns-nt)); % adapt window size

%% initialization
W=inf*ones(nt+1,ns+1); % cache matrix
W(1,1)=0;
if ~a
    W(1,:) = 0;
end

%% begin dynamic programming
for i=2:ns+1
    for j=max(i-w,1):min(i+w,nt+1)
        if j > 1
            oost = pdist2(s(i-1,:),t(j-1,:));    % Euclidean distance
%             oost=sum(abs(s(i-1,:)-t(j-1,:)));    % Absolut distance
            W(j,i)=(oost+min([W(j-1,i),W(j,i-1),W(j-1,i-1)]));
        end        
    end
end
% d=W(end,end);
