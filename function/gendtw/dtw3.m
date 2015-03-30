% Computes the dynamic time warping matrix of two multi-dimensional signals
% using dissimilarity cost matrix as the distance measure

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

function W=dtw3(s,t,a,w,D,KM,KT)
% Output:
% Input:
    % s: signal 1
    % t: signal 2
    % a: flag for alignment '1' or non-alignment (detection) '0'
    % w: window parameter
    %      if s(i) is matched with t(j) then |i-j|<=w
    % d: resulting distance
    % D: Dissimilarity matrix
    % KM: Subgesture cost for the model sequence
    % KT: Subgesture costs for the test sequence

if nargin<5
    D = [];
    if nargin<4
        w=Inf;
        if nargin<3
            a=true;
        elseif ~islogical(a)
            error('a parameter must be a logical value');
        end
    end
% elseif isempty(KM) || isempty(KT)
%     error('Matrices must be filled to compute the costs');
end

if size(s,2) > 1 && size(t,2) > 1
    ns=size(s,1);
    nt=size(t,1);
else
    ns=length(s);
    nt=length(t);
end
w=max(w, abs(ns-nt)); % adapt window size

%% initialization
W=zeros(nt+1,ns+1)+Inf; % cache matrix
W(1,1)=0;
if ~a
    W(1,:) = 0;
end

%% begin dynamic programming
for i=2:ns+1
    for j=max(i-w,1):min(i+w,nt+1)
%         W(i+1,j+1)=oost+min( [W(i,j+1), W(i+1,j), W(i,j)] );
%         oost=abs(s(i)-t(j));
        if j > 1
            if isempty(D)
                oost = pdist2(s(i-1,:),t(j-1,:));    % Euclidean distance
            else
%                 [~,pM] = min(KM(:,j-1));
%                 [~,pT] = min(KT(:,i-1));
                oost = D(i-1,j-1);    % Euclidean distance
%                 oost = D(pM,pT);
            end 
            W(j,i)=(oost+min([W(j-1,i),W(j,i-1),W(j-1,i-1)]));
        end        
    end
end
% d=W(end,end);
