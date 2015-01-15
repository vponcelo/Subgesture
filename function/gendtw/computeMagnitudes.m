function [segs,fr_fixed,params] = computeMagnitudes(X,params)
% Output:
%   N: Number of subsequences
%   segs: sequence segments
%   fr_fixed: mean length of subsequences
% Input:
%   X: training data sequences
%   params: struct with the required parameters

global COORDS;
global BASELINE;
global JOINTS;

if strcmp(COORDS,'world');
    ncoords = 3;
elseif strcmp(COORDS,'pixel')
    ncoords = 2;
end
Xgrad = zeros(size(X,1),ncoords*2+2);
Xhess = zeros(size(X,1),ncoords*2+2);
for i = 1:size(X,1)
    for coord = 0:ncoords-1
        [leftWrist,poslw] = ismember(7,JOINTS);
        if ~leftWrist
            error('computeMagnitudes:missJoint','Left wrist joint is not contained in the data');
        end
        [rightWrist,posrw] = ismember(11,JOINTS);
        if ~rightWrist
            error('computeMagnitudes:missJoint','Right wrist joint is not contained in the data');
        end
        if i == 1            
            Xgrad(1,coord+1) = X(1,poslw+length(JOINTS)*coord);            
            Xgrad(1,coord+1+ncoords) = X(1,posrw+length(JOINTS)*coord);
            Xhess(1,coord+1) = Xgrad(1,coord+1);
            Xhess(1,coord+1+ncoords) = Xgrad(1,ncoords+coord+1);
        else
            Xgrad(i,coord+1) = X(i,poslw+length(JOINTS)*coord)-X(i-1,poslw+length(JOINTS)*coord);
            Xgrad(i,coord+1+ncoords) = X(i,posrw+length(JOINTS)*coord)-X(i-1,posrw+length(JOINTS)*coord);
            Xhess(i,coord+1) = Xgrad(i,coord+1)-Xgrad(i-1,coord+1);
            Xhess(i,coord+1+ncoords) = Xgrad(i,ncoords+coord+1)-Xgrad(i-1,ncoords+coord+1);
        end
    end
end
Xgrad(:,end-1) = sqrt(sum(Xgrad(:,1:ncoords).^2')) + sqrt(sum(Xgrad(:,ncoords+1:end-2).^2'));
t_flow = min(Xgrad(:,end-1)) + (max(Xgrad(:,end-1))-min(Xgrad(:,end-1)))/(length(X)/params.thMinMov);   
% good values of t_flow for 3D: 0.0048 con params.thMinMov = 1
% good values of t_flow for 2D: 3.8962 con params.thMinMov = 1
% good values of t_flow for 2D: 2.9971 fixed con params.thMinMov = 1.3

Xgrad(:,end) = Xgrad(:,end-1) < t_flow;
Xhess(:,end-1) = sqrt(sum(Xhess(:,1:ncoords).^2')) == 0;
Xhess(:,end) = sqrt(sum(Xhess(:,ncoords+1:end-2).^2')) == 0;
fr_anchors =  Xgrad(:,end) | Xhess(:,end-1) | Xhess(:,end);

segs = []; c_ones = 0; ant_ones = false;
for i = 1:length(fr_anchors)
    if fr_anchors(i) == 0
        if i == 1 || i == length(fr_anchors)
            segs = [segs i];
        end
        if c_ones > 0 && ~ant_ones
            segs = [segs i-1];
            segs = [segs i];                      
        elseif ant_ones 
            segs = [segs i];
            ant_ones = false;
        end
        c_ones = 0;
    else
        c_ones = c_ones + 1;
        if c_ones == params.threshMov && ~ant_ones
            segs = [segs i];
            c_ones = 0;
            ant_ones = true;
        end
    end
end

if isempty(segs)
    error('computeMagnitudes:segNull','There is no initial segmentation');
end
                
if mod(length(segs),2) == 1
    if fr_anchors(1) == 1
        segs(1) = [];
    else
        segs(end) = [];
    end
end

fr_fixed = 0;
if ~strcmp(params.Baseline,BASELINE{1})    
    if strcmp(params.Baseline,BASELINE{3})
        params.N = length(segs)-1;
        lengths = zeros(1,params.N);
        for i = 1:length(segs)
            lengths(i) = segs(2*(i-1)+1)-segs(2*(i-1)+1)+1;
        end
        params.nmin = 1;
        params.nmax = max(lengths)+1;
        fr_fixed = round(mean(lengths));
        if ~fr_fixed
            error('computeMagnitudes:fixFrNull','There is no fixed subsequence length');
        end
        params.N = round(size(X,1)/fr_fixed);
    end
end

%% Show skeletons (only for 2D pixels)
if strcmp(COORDS,'pixel') && params.drawMovSkels
    drawAllSkels(X,segs,Xgrad(:,end),Xhess(:,end-1),Xhess(:,end))
end
