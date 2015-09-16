%% Global variables
global DATATYPE;        % datasets: 'random' 'chalearn2013' 'chalearn2014' 'madX' 'msr3dX/S' 'msract3d'
global VISUALIZE;       % flag indicating whether or not to visualize the data groups using kmeans clustering
global FRAMECOMPRESS;   % flag indicating whether or not to perform frame compression 
global KMEANSDTWv;      % possible versions of the k-means DTW algorithm
global DISTANCES;       % possible distances: 'dtw' 'euclidean'
global SAMPLING;        % sampling type
global S;               % Scores from the GA training
global Stest;           % Scores of the GA optimization in test
global CACHE;           % Cache for the GA
global BESTIND;         % Best individual at each generation
global PREDICTIONS;     % Predictions reached: [startFrame,endFrame,predLabel]
global BASELINE;        % Baseline to be used in the GA
global COORDS;          % type of coordinates
global STATE;           % save state of previous GA run
global OPTIONS;         % Options of the GA run
global PERCENTDATA;     % percent of learning data to consider
global NORMTYPE;        % Normalization type 'neck' 'none'
global MEDIANTYPE;      % Type of median models 
global JOINTS;          % Selected joints
global NAT;             % Type of Descriptor

DATATYPE = 'chalearn2013';
NORMTYPE = 'neck';
COORDS = 'pixel';
VISUALIZE = false;
FRAMECOMPRESS = true;
KMEANSDTWv = {'v1_0','v1_1','v1_2','v1_3','v2_0_0','v2_0_1','v2_1_0',...
    'v2_1_1','v2_2_0','v2_2_1','v2_3_0','v2_3_1'};
DISTANCES = {'dtwCost','euclidean'};
SAMPLING = {'random','labels','segments'};
S = []; Stest = [0 0];
BASELINE = {'random','deriv','fixed'};
STATE = [];
BESTIND.model = []; BESTIND.state = [];
OPTIONS = [];
PERCENTDATA = 100;
MEDIANTYPE = {'direct','modelSM1','modelSM2','allSM1','allSM2','DCSR'};
JOINTS = [4 6 7 8 10 11 12];%1:20;%  %[8,12] hands
NAT = 0;

%% parameters data structures
nrsamples = 100;        % number of random samples
nseqs = 0;              % number of batches for training and validation, respectively. '0': All ; [nTrain nVal]
nframesSeg = 0;         % Fixed number of frames for subgesturing. '0' means no fixed subgesturing
sampling = SAMPLING{2}; % sampling type: 'random' 'label' 'segments'
noise = true;          % flag indicating whether or not consider noise (iddle gesture) for the test sequence
secsBatch = 60;         % reference seconds for the test sequence
nSampGest = 0;        % Number of samples per gesture for the test sequence

%% classification type parameter
params.classification = true;       % flag that indicates to perform global category classification
params.accuracyglobal = true;        % global accuracy or weighted accuracy (for imbalanced data sets).
params.darwin = true;                % flag that indicates to perform videoDarwin-based classification
params.svm = false;                  % flag that indicates to perform SVM-based classification
params.tc = false;                   % flag that indicates whether to perform temporal clustering

%% parameters hmm
params.phmm.folds = 1;                 % k for k-fold Cross Validation
params.phmm.states = 3;                % number of hidden states for the HMM
params.phmm.it = 1000;                  % number of Iterations of the HMM
params.phmm.clustType = 'none';        % clustering method: 'none' 'kmlsample' 'kmeans' 'haca'
params.phmm.kD = 300;                   % number of clusters for discretizing
params.phmm.cIters = 100;              % number of iterations for discretizing
params.phmm.varType = 'discrete';  % type of variable for the HMM: 'gauss' 'mixgausstied' 'discrete' 
params.phmm.hmm = true;            % flag that indicates to train with hmm training 
params.phmm.pmtk = true;           % flag that indicates to use the pmtk3 library implementation

%% parameters genetic temporal clustering
params.version = ...
    char(KMEANSDTWv{5});    % versions of the k-means DTW algorithm to execute';
params.dist = DISTANCES{1}; % distance metric
params.k0 = 3;              % initial data clusters for subgesturing (< 0: 'no clustering')
params.nmin = 5;            % minimum subsequence width
params.nmax = 25;           % maximum subsequence width
params.N = 500;             % Number of segments to split the learning sequence
params.N0 = 8;              % Number of segments to split the model sequences
params.nThreshs = 20;       % Number of thresholds for testing (tunned to 20 and 22 for max models and median models, respectively)
params.D = [];              % Dissimilarity matrix
params.bestThs = [];        % Thresholds learnt on training
params.vectorized = 'off';   % vectorize the GA
params.population = 20;     % population of the GA
params.generations = 1006;  % number of generations of the GA
params.Baseline = ...
    BASELINE{2};            % Baseline for the GA
params.threshMov = 3;       % maximum number of low movement frames
params.thMinMov = 1.3;      % Minimum portion of movement
params.drawMovSkels = false;% flag indicating whether to draw skels
params.probMut = 0.6;         % probability of performing standard or specific GA mutation '1':= only standard
params.probCross = 0.6;     % probability of performing standard or specific GA crossover '1':= only standard
params.scale = 0.5;         % scale parameter for Gaussian mutation
params.shrink = 0.75;       % shrink parameter for Gaussian mutation
params.probSeg = 0.2;       % probability of eliminate/change a segment
params.maxWlen = 1000;      % maximum DTW cost matrix length to detect the start-end
params.mCostType = 'mean';  % 'mean' subgesture costs for the models
params.mType = MEDIANTYPE{3};  % Type of median models to consider
params.usemax_l = true;        % use the median or the max-length gesture as reference
params.resize = true;          % Use resizing instead of mean DTW alignment
params.gmm = false;            % Use gmm instead of other non-probabilistic representations
params.pdtw = false;           % flag for indicating the use of gmms in feature modeling
params.score2optim = 'p';        % Score to optimize --> Overlap: 'o', Precision: 'p', Recall: 'r', Accuracy/F1-Score(spotting): 'a'
params.minOverlap = 0.5;        % Minimum overlap to detect the label
params.sw = 0;              % sliding window (frame seq length): '0' means the whole sequence
params.k = 0;               % current k to evaluate for the K-Nearest Neighbour DTW models
CACHE.pos = int32(1);       % Index positions
% GENRESULTS.P = cell(1,params.generations);
% GENRESULTS.eval = cell(1,params.generations);
PREDICTIONS = [];