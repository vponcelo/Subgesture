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
global NAT;             % Type of Descriptor
global MSM;             % Type of Median Subgesture Models to use: 'fixed' or 'evolutive'

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
MEDIANTYPE = {'direct','modelMSM1','modelMSM2','allMSM1','allMSM2','KNN','DCSR'};
JOINTS = [4 6 7 8 10 11 12];%1:20;%  %[8,12] hands
NAT = 0;
MSM = {'none','fix','evoSegs'};

%% parameters data structures
nrsamples = 100;        % number of random samples
nseqs = 0;              % number of batches for training and validation, respectively. '0': All ; [nTrain nVal]
nframesSeg = 0;         % Fixed number of frames for subgesturing. '0' means no fixed subgesturing
sampling = SAMPLING{2}; % sampling type: 'random' 'label' 'segments'
noise = true;          % flag indicating whether or not consider noise (iddle gesture) for the test sequence
secsBatch = 60;         % reference seconds for the test sequence
nSampGest = 0;          % Number of samples per gesture for the test sequence   

%% parameters hmm
params.phmm.folds = 1;                 % k for k-fold Cross Validation
params.phmm.states = 3;                % number of hidden states for the HMM
params.phmm.it = 1000;                  % number of Iterations of the HMM
params.phmm.clustType = 'none';        % clustering method: 'none' 'kmlsample' 'kmeans' 'haca'
params.phmm.kD = 50;                   % number of clusters for discretizing
params.phmm.cIters = 100;              % number of iterations for discretizing
params.phmm.varType = 'discrete';  % type of variable for the HMM: 'gauss' 'mixgausstied' 'discrete' 
params.phmm.hmm = false;            % flag that indicates to train with hmm training 

%% parameters genetic temporal clustering
params.version = ...
    char(KMEANSDTWv{5});    % versions of the k-means DTW algorithm to execute';
params.dist = DISTANCES{1}; % distance metric
params.k0 = 3;              % initial data clusters for subgesturing
params.nmin = 5;            % minimum subsequence width
params.nmax = 25;           % maximum subsequence width
params.N = 500;             % Number of segments to split the learning sequence
params.N0 = 8;              % Number of segments to split the model sequences
params.nThreshs = 20;       % Number of thresholds for testing (tunned to 20 and 22 for max models and median models, respectively)
params.D = [];              % Dissimilarity matrix
params.bestThs = [];        % Thresholds learnt on training
params.vectorized = 'off';   % vectorize the GA
params.population = 10;     % population of the GA
params.generations = 1000;  % number of generations of the GA
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
params.msmType = MSM{1};    % Type of Median Subgesture Models in the evolutive process 
params.mType = MEDIANTYPE{3};  % Type of median models to consider
params.usemax_l = true;        % use the median or the max-length gesture as reference
params.resize = true;          % Use resizing instead of mean DTW alignment
params.gmm = false;            % Use gmm instead of other non-probabilistic representations
params.sw = 5000;              % sliding window (frame seq length): '0' means the whole sequence
if strcmp(params.mType,'KNN')
	params.k = 3;
else
	params.k = 0;           % current k to evaluate for the K-Nearest Neighbour DTW models
end
CACHE.pos = int32(1);       % Index positions
% GENRESULTS.P = cell(1,params.generations);
% GENRESULTS.eval = cell(1,params.generations);
PREDICTIONS = [];