% 
% (c) 2019 Naoki Masuyama
% 
% These are the codes of Topological CIM-based Adaptive Resonance Theory (TCA)
% proposed in "N. Masuyama, C. K. Loo, H. Ishibuchi, N. Kubota, Y. Nojima, 
% and Y. Liu, "Topological clustering via adaptive resonance theory with information 
% theoretic learning," IEEE Access, vol. 7, no. 1, pp. 76920-76936, December 2019."
% 
% Please contact "masuyama@cs.osakafu-u.ac.jp" if you have any problems.
%     


MIter = 5; % Number of iterations

% Load Data
load 2D_2Gauss_1Uniform_Noise10


% Parameters of TCA ===================================================
TCAnet.numNodes    = 0;   % Number of clusters
TCAnet.weight      = [];  % Mean of cluster
TCAnet.CountNode = [];    % Counter for each node
TCAnet.edge = zeros(2,2); % Initial connections (edges) matrix
TCAnet.NewEdgedNode = []; % Node which creates new edge.
TCAnet.ErrCIM = [];       % CIM between nodes
TCAnet.adaptiveSig = [];
TCAnet.diligence   = [];
TCAnet.age = zeros(2,2);

TCAnet.Lambda = 200;      % Interval for Node deletion and topology construction
TCAnet.initSig = 0.12;    % Initial Kernel Bandwidth for CIM
% ====================================================================


for nitr = 1:MIter
    fprintf('Iterations: %d/%d\n',nitr,MIter);
    
    % Randamize data
    ran = randperm(size(DATA,1));
    DATA = DATA(ran,:);
    
    TCAnet = TCA(DATA, TCAnet);
    plotTCA(DATA, TCAnet);
    
end


