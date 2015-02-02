%% Global variables
global DATATYPE;        % dataset: 'chalearn' 'random'
global VISUALIZE;       % flag indicating whether or not to visualize the data groups using kmeans clustering
global FRAMECOMPRESS;   % flag indicating whether or not to perform frame compression 
%global ALLLEARNING;     
global KMEANSDTWv;      % possible versions of the k-means DTW algorithm
global DISTANCES;       % possible distances: 'dtw' 'euclidean'
global SAMPLING;        % sampling type
global S;               % Scores from the GA training
global CACHE;           % Cache for the GA
% global GENRESULTS;      % Population at each generation and results
global PREDICTIONS;     % Predictions reached: [startFrame,endFrame,predLabel]
global BASELINE;        % Baseline to be used in the GA
global COORDS;          % type of coordinates
global STATE;           % save state of previous GA run
global OPTIONS;         % Options of the GA run
global PERCENTDATA;     % percent of learning data to consider
global NORMTYPE;        % Normalization type 'neck' 'xyzangles'
global MEDIANTYPE;      % Type of median models 
global JOINTS;          % Selected joints
global NAT;             % Descriptor

DATATYPE = 'chalearn2014';
NORMTYPE = 'neck';
COORDS = 'pixel';
VISUALIZE = false;
FRAMECOMPRESS = true;
%ALLLEARNING = true;
KMEANSDTWv = {'v1_0','v1_1','v1_2','v1_3','v2_0_0','v2_0_1','v2_1_0',...
    'v2_1_1','v2_2_0','v2_2_1','v2_3_0','v2_3_1'};
DISTANCES = {'dtwCost','euclidean'};
SAMPLING = {'random','labels','segments'};
S = [];
BASELINE = {'random','deriv','fixed'};
STATE = [];
OPTIONS = [];
PERCENTDATA = 100;
MEDIANTYPE = {'direct','directMSM1','directMSM2','KNN','DCSR'};
JOINTS = [4 6 7 8 10 11 12];%1:20;%  %[8,12] hands
NAT = 0;

%% parameters data structures
nrsamples = 100;        % number of random samples
nseqs = 0;              % number of batches for training and validation, respectively. '0': All ; [nTrain nVal]
nframesSeg = 0;         % Fixed number of frames for subgesturing. '0' means no fixed subgesturing
sampling = SAMPLING{2}; % sampling type: 'random' 'label' 'segments'
noise = false;          % flag indicating whether or not consider noise (iddle gesture) for the test sequence
secsBatch = 60;         % reference seconds for the test sequence
nSampGest = 0;          % Number of samples per gesture for the test sequence   

% classification
folds = 1;              % k for k-fold Cross Validation
clustType = 'kmeans';   % clustering method: 'kmlsample' 'kmeans' 'haca'
numClusters = 100;       % number of clusters for discretizing
numIterations = 100;    % number of iterations for discretizing
hmmIters = 50000;       % number of Iterations of the HMM

% genetic temporal clustering parameters
params.version = ...
    char(KMEANSDTWv{5});        % versions of the k-means DTW algorithm to execute';
params.dist = DISTANCES{1}; % distance metric
params.k0 = 3;              % initial data clusters for subgesturing
params.nmin = 5;            % minimum subsequence width
params.nmax = 25;           % maximum subsequence width
params.N = 500;             % Number of segments to split the learning sequence
params.N0 = 5;             % Number of segments to split the model sequences
params.nThreshs = 100;      % Number of thresholds for testing
params.D = [];              % Dissimilarity matrix
params.bestThs = [];        % Thresholds learnt on training
params.vectorized = 'on';   % vectorize the GA
params.population = 5;      % population of the GA
params.generations = 2;  % number of generations of the GA
params.Baseline = ...
    BASELINE{2};            % Baseline for the GA
params.threshMov = 3;       % maximum number of low movement frames
params.thMinMov = 1.3;      % Minimum portion of movement
params.drawMovSkels = false;% flag indicating whether to draw skels
params.thMutCross = 0.6;    % threshold of GA mutation and crossover
params.scale = 0.5;         % scale parameter for Gaussian mutation
params.shrink = 0.75;       % shrink parameter for Gaussian mutation
params.probSeg = 0.2;       % probability of eliminate/change a segment
params.maxWlen = 1000;      % maximum DTW cost matrix length to detect the start-end
params.msm = true;         % use Median Subgesture Models in the evolutive process instead of Median Models
params.mType = MEDIANTYPE{3};  % Type of median models to consider
if strcmp(params.mType,'KNN')
	params.k = 3;
else
	params.k = 0;           % current k to evaluate for the K-Nearest Neighbour DTW models
end
CACHE.ind = int32(zeros(... % Cache with the populations
    params.population*100,params.N*2+1,'int32'));
CACHE.pos = int32(1);       % Index positions
% GENRESULTS.P = cell(1,params.generations);
% GENRESULTS.eval = cell(1,params.generations);
PREDICTIONS = [];