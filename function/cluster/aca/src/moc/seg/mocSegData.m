function wsData = mocSegData(wsSrc, varargin)
% Prepare motion capture data for temporal segmentation.
%
% Input
%   wsSrc    -  mocap source
%
% Output
%   wsData
%     sR     -  temporal reduction
%     Dof    -  Degree of Freedom, dim x nF0
%     X      -  quaternion (continuous) feature, dim x nF
%     X0     -  quaternion (continuous) feature without temporal reduction, dim x nF0
%     XD     -  1-D (discrete) feature, 1 x nF
%     XD0    -  1-D (discrete) feature without temporal reduction, 1 x nF0
%     segT   -  ground truth segmentation
%     segT0  -  ground truth segmentation without temporal reduction
%     cnames -  class names
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 12-29-2008
%   modify   -  Feng Zhou (zhfe99@gmail.com), 01-04-2010

m = length(varargin);

% source
[src, paraF] = stFld(wsSrc, 'src', 'paraF');

if m,
    % feature
    X0 = varargin{1};
    segT0 = varargin{2};
    cnames = varargin{3};
else
    prom('t', 'new moc seg data\n');
    
    % Dof
    wsDof = mocDof(src);

    % feature
    wsFeat = mocFeat(src, paraF, wsDof);
    X0 = wsFeat.X;
end

% temporal reduction
[sR, X, XD, XD0] = tempoReduce(X0, wsSrc.para);

if ~m,    
    % ground truth
    MOC = mocHuman;
    segT0 = MOC.(src.name).seg;
    cnames = MOC.(src.name).cnames;
end

% store
wsData.sR = sR;
wsData.X = X;
wsData.XD = XD;
wsData.X0 = X0;
wsData.XD0 = XD0;
wsData.segT0 = segT0;
wsData.segT = segNewS(segT0, sR, 'type', 'shrink');
wsData.cnames = cnames;
