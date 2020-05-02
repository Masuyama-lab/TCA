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
function net = TCA(DATA, net)


numNodes = net.numNodes;       % Number of nodes
weight = net.weight;           % Mean of node
CountNode = net.CountNode;     % Counter for each node
initSig = net.initSig;   % Initial Kernel Bandwidth for CIM
adaptiveSig = net.adaptiveSig; % cimSig in each node
diligence = net.diligence;

% Parameters for Topology
edge = net.edge;
NewEdgedNode = net.NewEdgedNode;
ErrCIM = net.ErrCIM;
Lambda = net.Lambda;


    
for sampleNum = 1:size(DATA,1)
    
    
    % Current data sample
    input = DATA(sampleNum,:);
    
    % Similarity measurement between an input and nodes.
    nodeCIM = CIM(input, weight, mean(adaptiveSig));
    [stateCIM, orderCIM] = sort(nodeCIM, 'ascend');
    
    % Number of nodes close to input in the range of mean(adaptiveCIM)
    NstateCIM = sum( stateCIM <= mean(adaptiveSig) );
    
    % for resonance process
    resonance = false;
    currentSortedIndex = 1;
    
    
    
    if NstateCIM == 0 % No clusters around input
        
        % If there is no node in the space.
        if isempty(weight) == 1
            adaptiveSig(1, 1) = initSig;
        end
        
        % Add Node
        numNodes                 = numNodes + 1;
        weight(numNodes,:)       = input;
        CountNode(1, numNodes)   = 1;
        adaptiveSig(1, numNodes) = mean(adaptiveSig);
        diligence(1, numNodes)   = 1;
        
        NewEdgedNode(1, numNodes) = 0;
        ErrCIM(1, numNodes)       = 1;
        edge(numNodes, :)         = 0;
        edge(:, numNodes)         = 0;
        
    elseif NstateCIM >= 1
        
        while ~resonance
            
            % Find a winner node
            s1 = orderCIM(currentSortedIndex);
            
            % Update the winner node
            bestWeight = (CountNode(1,s1)/(CountNode(1,s1)+1))*weight(s1,:) + (1/(CountNode(1,s1)+1))*input;
            
            % Calculate CIM between the winner node and an input for Vigilance Test
            bestCIM = CIM(input, bestWeight, adaptiveSig(1, s1));
                
            % Vigilance Test
            if bestCIM <= mean(adaptiveSig)
                % Match Success   
                % Update Parameters
                weight(s1,:)     = bestWeight;
                CountNode(1, s1) = CountNode(1, s1) + 1;
                diligence(1, s1) = diligence(1, s1) + 1;
                
                resonance = true;
                  
                % Topology Construction
                % Calculate CIM based on t-th and (t+1)-th node position as an Error state.
                Ecim = CIM(weight(s1,:), bestWeight, adaptiveSig(1, s1)); % If CountCluster is large, Ecim goes to zero.
                
                % If an Error become small comparing with previous state, update ErrCIM.
                ErrCIM(1, s1) = min(ErrCIM(1, s1), Ecim);
                
                % Create an edge between s1 and s2 nodes.
                if NstateCIM >= 2
                    s2 = orderCIM(currentSortedIndex + 1); % Second Matching node Index
                    if CIM(input, weight(s2,:), adaptiveSig(1, s2)) <= mean(adaptiveSig) % Vigilance Test for s2 node.
                        edge(s1,s2) = 1;
                        edge(s2,s1) = 1;
                        NewEdgedNode(1,s1) = 1;
                    end
                end
                
            else
                % Match Fail
                if(currentSortedIndex == numNodes)  % Reached to maximum number of generated nodes.
                    % Add Node
                    [idx,~] = knnsearch(weight,input,'k',3); % Added node should have a similar state to neighbors.
                    numNodes                 = numNodes + 1;
                    CountNode(1, numNodes)   = mean(CountNode(1, idx));
                    adaptiveSig(1, numNodes) = mean(adaptiveSig(1, idx));
                    weight(numNodes,:)       = input;
                    diligence(1, numNodes)   = 1;
                    diligence(1, :)          = zeros(size(numNodes)); % Once added a new node, reset deligence.
                    
                    NewEdgedNode(1, numNodes) = 0;
                    ErrCIM(1, numNodes)       = 1;
                    edge(numNodes, :)         = 0;
                    edge(:, numNodes)         = 0;
                    
                    resonance = true;
                else
                    currentSortedIndex = currentSortedIndex + 1;    % Search another cluster orderd by sortedNodes.
                end

            end % end vigilance

        end % end resonance
        
    end % NstateCIM test
    
    
    % Topology Reconstruction
    if mod(sampleNum, Lambda) == 0
        
        % -----------------------------------------------------------------
        % Delete Node based on number of neighbors
        nNeighbor = sum(edge);
        deleteNodeEdge = (nNeighbor == 0);
        
        if sum(deleteNodeEdge) ~= size(weight,1)
            % Delete process
            edge(deleteNodeEdge, :) = [];
            edge(:, deleteNodeEdge) = [];
            weight(deleteNodeEdge, :) = [];
            numNodes = numNodes - sum(deleteNodeEdge);
            CountNode(:, deleteNodeEdge) = [];
            NewEdgedNode(:, deleteNodeEdge) = [];
            ErrCIM(:, deleteNodeEdge) = [];
            adaptiveSig(:, deleteNodeEdge) = [];
            diligence(:, deleteNodeEdge) = [];
        end
        
        % -----------------------------------------------------------------
        % Delete Node based on ErrCIM
        [stateEC, posEC] = sort(ErrCIM, 'ascend');
        highEC = ( stateEC > mean(adaptiveSig)*2);
        deleteNodeEC = posEC(highEC);
        
        if size(deleteNodeEC,2) ~= size(weight,1)
            % Delete process
            edge(deleteNodeEC, :) = [];
            edge(:, deleteNodeEC) = [];
            weight(deleteNodeEC, :) = [];
            numNodes = numNodes - size(deleteNodeEC,2);
            CountNode(:, deleteNodeEC) = [];
            NewEdgedNode(:, deleteNodeEC) = [];
            ErrCIM(:, deleteNodeEC) = [];
            adaptiveSig(:, deleteNodeEC) = [];
            diligence(:, deleteNodeEC) = [];
        end
        
        % -----------------------------------------------------------------
        % Delete Intersections of edge
        [weight, edge, NewEdgedNode] = DeleteIntersection(weight, edge, NewEdgedNode, mean(adaptiveSig));
        
    end % end topology reconstruction
    
    
    
    % Updates cimSig if there is no node inserted for a certain period.
    if min(diligence) == 0
        cnd = 1;
    else
        cnd = min(diligence);
    end
    
    if mean(diligence) > cnd*10 && numNodes > 1
        % Cluster Labeling based on edge
        connection = graph(edge ~= 0);
        LebelCluster = conncomp(connection);
        for k = 1:max(LebelCluster)
            idxLC = find(LebelCluster == k);
            selectedDili = diligence(idxLC);

            % Check uniformity of data distribution
            % DATA has an uniform distribution -> uniformDegree = 0.5
            uniformDegree = median(selectedDili)/(median(selectedDili) + mean(selectedDili));
            
            % If DATA does not have an uniform distribution, cimSig is updated.
            if (0.5 - uniformDegree) > 0.05
                [~, upper] = find(selectedDili >= mean(selectedDili)); % for Nodes in the dense region
                [~, lower] = find(selectedDili < mean(selectedDili));  % for Nodes in the sparse region
                for i = upper
                    rate = abs(uniformDegree)*0.1;
                    adaptiveSig(1, idxLC(i)) = adaptiveSig(1, idxLC(i)) - rate * mean(adaptiveSig(idxLC(upper)));
                    if adaptiveSig(1, idxLC(i)) <= initSig*0.5
                        adaptiveSig(1, idxLC(i)) = initSig*0.5;
                    end
                end
                for j = lower
                    rate = abs(uniformDegree)*0.1;
                    adaptiveSig(1, idxLC(j)) = adaptiveSig(1, idxLC(j)) + rate * mean(adaptiveSig(idxLC(lower)));
                    if adaptiveSig(1, idxLC(j)) >= initSig*1.5
                        adaptiveSig(1, idxLC(j)) = initSig*1.5;
                    end
                end
            elseif (0.5 - uniformDegree) <= -0.05
                [~, upper] = find(selectedDili >= mean(selectedDili)); % for Nodes in the dense region
                [~, lower] = find(selectedDili < mean(selectedDili));  % for Nodes in the sparse region
                for i = upper
                    rate = abs(uniformDegree)*0.1;
                    adaptiveSig(1, idxLC(i)) = adaptiveSig(1, idxLC(i)) + rate * mean(adaptiveSig(idxLC(upper)));
                    if adaptiveSig(1, idxLC(i)) >= initSig*1.5
                        adaptiveSig(1, idxLC(i)) = initSig*1.5;
                    end
                end
                for j = lower
                    rate = abs(uniformDegree)*0.1;
                    adaptiveSig(1, idxLC(j)) = adaptiveSig(1, idxLC(j)) - rate * mean(adaptiveSig(idxLC(lower)));
                    if adaptiveSig(1, idxLC(j)) <= initSig*0.5
                        adaptiveSig(1, idxLC(j)) = initSig*0.5;
                    end
                end
            end
            
        end

        % Reset diligence
        diligence(1, :) = zeros(size(numNodes));

    end

end % end numSample


% -----------------------------------------------------------------
% Delete Node based on number of neighbors
nNeighbor = sum(edge);
deleteNodeEdge = (nNeighbor == 0);

if sum(deleteNodeEdge) ~= size(weight,1)
    % Delete process
    edge(deleteNodeEdge, :) = [];
    edge(:, deleteNodeEdge) = [];
    weight(deleteNodeEdge, :) = [];
    numNodes = numNodes - sum(deleteNodeEdge);
    CountNode(:, deleteNodeEdge) = [];
    NewEdgedNode(:, deleteNodeEdge) = [];
    ErrCIM(:, deleteNodeEdge) = [];
    adaptiveSig(:, deleteNodeEdge) = [];
    diligence(:, deleteNodeEdge) = [];
end

% -----------------------------------------------------------------
% Delete Node based on ErrCIM
[stateEC, posEC] = sort(ErrCIM, 'ascend');
highEC = ( stateEC > mean(adaptiveSig)*2);
deleteNodeEC = posEC(highEC);

if size(deleteNodeEC,2) ~= size(weight,1)
    % Delete process
    edge(deleteNodeEC, :) = [];
    edge(:, deleteNodeEC) = [];
    weight(deleteNodeEC, :) = [];
    numNodes = numNodes - size(deleteNodeEC,2);
    CountNode(:, deleteNodeEC) = [];
    NewEdgedNode(:, deleteNodeEC) = [];
    ErrCIM(:, deleteNodeEC) = [];
    adaptiveSig(:, deleteNodeEC) = [];
    diligence(:, deleteNodeEC) = [];
end

% -----------------------------------------------------------------
% Delete Intersections of edge
[weight, edge, NewEdgedNode] = DeleteIntersection(weight, edge, NewEdgedNode, mean(adaptiveSig));


% Cluster Labeling based on edge
connection = graph(edge ~= 0);
LebelCluster = conncomp(connection);



net.numNodes = numNodes;      % Number of nodes
net.weight = weight;          % Mean of nodes
net.CountNode = CountNode;    % Counter for each node
net.adaptiveSig = adaptiveSig;
net.diligence = diligence;

net.LebelCluster = LebelCluster;
net.edge = edge;
net.NewEdgedNode = NewEdgedNode;
net.ErrCIM = ErrCIM;

net.Lambda = Lambda;

end



% Correntropy induced Metric (Gaussian Kernel based)
function cim = CIM(X,Y,sig)
% X : 1 x n
% Y : m x n
[n, att] = size(Y);
g_Kernel = zeros(n, att);

for i = 1:att
    g_Kernel(:,i) = GaussKernel(X(i)-Y(:,i), sig);
end

ret0 = GaussKernel(0, sig);
ret1 = mean(g_Kernel, 2);

cim = sqrt(ret0 - ret1)';
end


function g_kernel = GaussKernel(sub, sig)
g_kernel = exp(-sub.^2/(2*sig^2));
% g_kernel = 1/(sqrt(2*pi)*sig) * exp(-sub.^2/(2*sig^2));
end



% Delete intersections of edge
function [weight, edge, NewEdgedNode] = DeleteIntersection(weight, edge, NewEdgedNode, sigma)

% for d = 1:size(weight,1); % Search all nodes
for d = find(NewEdgedNode == 1) % Search only new edged nodes
    
    node1 = find(edge(d,:)); % Neighbors of d-th node
    if size(node1,1) >= 1
       posX1 = weight(d,:); % position of d-th node
        for m = 1:size(node1,2) % Search all neighbors of d-th nodes
            posY1 = weight(node1(m),:); % position of m-th neighbor node of d-th node
            for h = 1:size(node1,2)
                target2 = node1(h);
                node2 = find(edge(target2,:)); % Neighbors of m-th node
                posX2 = weight(target2,:); % position of h-th neighbor node of m-th node
                for k = 1:size(node2,2)
                    posY2 = weight(node2(k),:); % position of k-th neighbor node of h-th node
                    isConvex = findIntersection(posX1, posY1, posX2, posY2); % find intersections
                    if isConvex == 1 % If intersection is exist, delete edge which has larger CIM.
                        cim1 = CIM(weight(d,:), weight(node1(m),:), sigma);
                        cim2 = CIM(weight(target2,:), weight(node2(k),:), sigma);
                        if cim2 >= cim1
                            edge(target2, node2(k)) = 0;
                            edge(node2(k), target2) = 0;
                        else
                            edge(d, node1(m)) = 0;
                            edge(node1(m), d) = 0;
                        end
                    end % end isConvex
                end % end k
            end % end h
        end % end m  
    end

end % end d

NewEdgedNode = zeros(size(NewEdgedNode));

end

% Check intersection of edges
function [isConvex] = findIntersection(A, B, C, D)

F1  = B(:,1)-D(:,1);
F2  = B(:,2)-D(:,2);
M11 = B(:,1)-A(:,1);
M21 = B(:,2)-A(:,2);
M12 = C(:,1)-D(:,1);
M22 = C(:,2)-D(:,2);
deter = M11.*M22 - M12.*M21;
lambda = -(F2.*M12-F1.*M22)./deter;
gamma = (F2.*M11-F1.*M21)./deter;

% E = (lambda*[1 1]).*A + ((1-lambda)*[1 1]).*B;
% isConvex = (0 <= lambda & lambda <= 1)  & (0 <= gamma & gamma <= 1);

isConvex = (0 < lambda & lambda < 1)  & (0 < gamma & gamma < 1) ;
isConvex = isConvex';

end




