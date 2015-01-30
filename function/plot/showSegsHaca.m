function showSegmentation(segType)
% show the segmentations results for visual understandability. 

% input
%   varargin: set of variables when HACA is running

%global FRAMECOMPRESS;

%if FRAMECOMPRESS,
        load(strcat('haca',segType));
%    else
%        load(strcat('haca',segType));
%end
    

% if length(varargin) == 6,
%     K = varargin{1};
%     X = varargin{2};
%     segT = varargin{3};
%     segAca = varargin{4};
%     segHaca = varargin{5};
%     segGmm = varargin{6};
% else
%     
% end

%% plot
showM(K, 'fig', [1 1 2 1]);
title('Kernel matrix (K)');
showSeq(X, 'fig', [1 1 2 2]);
title('feature in 2-D space');
showSegBar(segT,   'fig', [2 4 1 1], 'mkSiz', 0, 'lnWid', 1);
showSegBar(segAca, 'fig', [2 4 1 2], 'mkSiz', 0, 'lnWid', 1);
title(sprintf('aca accuracy %.2f', segAca.acc));
showSegBar(segHaca(end), 'fig', [2 4 1 3], 'mkSiz', 0, 'lnWid', 1);
title(sprintf('haca accuracy %.2f', segHaca(end).acc));
showSegBar(segGmm, 'fig', [2 4 1 4], 'mkSiz', 0, 'lnWid', 1);
title(sprintf('gmm accuracy %.2f', segGmm.acc));