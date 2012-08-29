function [S, N, VI, C] = stability_new(G, T, varargin)
%STABILITY    Graph partitioning optimizing stability with the Louvain
%             algorithm
%
%   [S, N, VI, C] = STABILITY(G, T) finds the optimal partitions of the 
%   graph G by optimizing the stability at each Markov time in vector T. G 
%   can either be the list of the edges in the graph (in the form [node i, 
%   node j, weight of link i-j; node k, node l, weight of link k-l;...] if 
%   the graph is weighted, or [node i, node j; node k, node l;...] if the 
%   graph is unweighted) or the adjacendy matrix of the graph. S, N, VI and 
%   C contain respectively the stability, the number of communities, the 
%   variation of information, and the optimal partition for each 
%   Markov time contained in T. If T is not specified, the modularity
%   (equivalent to stability for T=1) is calculated. Ideally, Markov time
%   should be sampled exponentially (e.g.: T = 10.^[-2:0.01:2]).
%
%   [S, N, VI, C] = STABILITY(G, T,'PARAM',VALUE) accepts one or more
%   comma-separated parameter name/value pairs. For a list of parameters 
%   and values, see "Parameter Options."
%
%    
%   Parameter Options:
% 
%        Parameter      Value                                 Default
%        ---------      -----                                 -------
%        L              Number of optimisations of the          100
%                       Louvain algorithm to be done at 
%                       each Markov time.
%
%        M              The top M partitions among the L        100
%                       given at each Markov time by the
%                       L louvain optimisations will be
%                       used to compute the variation of
%                       information.    
%
%        laplacian      Allows to choose which type of     'normalised'
%                       laplacian should be used to 
%                       calculate the stability. It can
%                       either be 'combinatorial', or 
%                       'normalised'.
%
%        directed       activate stability for directed         none
%                       graphs. Note that transition matrices
%                       are defined for left multiplications 
%                       here, i.e. A_ij is the link from i to j.
%
%
%        noVI           Disables the calculation of the         none
%                       robustness of the partitions.
%                       Disabling this can significantly 
%                       speed up the calculations.
%
%        out            Enables saving step by step the         ''
%                       partitions found at each Markov 
%                       time in a file located in a 
%                       folder named 'Partitions', as 
%                       well as the outputs in a text 
%                       file output_'...'.stdout matlab 
%                       file output_'...'.mat
%
%       linearised      Enables the calculation of the          none
%                       linearised stability instead of the 
%                       full stability (default).
%
%       nocheck         Disables the checks for the             none
%                       encoding of the graph. This can 
%                       save computational time but can 
%                       also lead to serious errors if 
%                       the graph has not been properly
%                       encoded.
%
%       prec            Precision: defines a threshold for     10e-9  
%                       the range of weights allowed in  
%                       the laplacian exponential matrix 
%                       of the full stability.
%
%       plot            Plots the plots the results of          none
%                       the stability, number of 
%                       communities and variation of 
%                       information as a function of the
%                       Markov time.           
%
%       v               Verbose mode                            none
%
%       p               Parallel mode                           none
%
%       t               Output as text files:                   none
%                       Enables saving step by step the         
%                       partitions found at each Markov 
%                       time in a text file located in a 
%                       folder named 'Partitions', as 
%                       well as the outputs in a text 
%                       file output_'...'.stdout matlab 
%                       The option 'out' must be on.
%
%


% Unparsed default parameters
Graph = [];                                     % List of edges of the graph to be partitioned
Time = 1;                                       % Markov times at which the graph should be partitioned
flag_matlabpool = false;                        % for opening/closing workpool for parallel computation



%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%$                                          $%
%$          Arguments parsing               $%
%$                                          $%
%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%

[StabilityFunction, OutputFile, prefix, Full, Sanity, plotStability, verbose, TextOutput, PARAMS] = parseinput(length(varargin),varargin);



% Argument 1: Graph

if nargin > 0
    if size(G,1) == size(G,2) && ~issparse(G)
            G=sparse(G);
    end
    % Check if the graph is correctly encoded
    if Sanity
        G=check(G, verbose, PARAMS);
    end
    % If the full stability is to be computed, Graph should be the
    % adjacency matrix.
    if Full
        if size(G,1) ~= size(G,2)
            if size(G,2)==3
                Graph=sparse(G(:,1)+1,G(:,2)+1,G(:,3));
            elseif size(G,2)==2
                Graph=sparse(G(:,1)+1,G(:,2)+1,ones(length(G(:,1)),1));
            else
                error('Wrong size for G: G should be a graph saved either as a list of edges (size(G)=[N,3] if weighted, size(G)=[N,2] if unweighted) or as an adjacency matrix (size(G)=[N,N])');
            end
        else
            Graph = sparse(G);
        end
        PARAMS.NbNodes = size(Graph,2);
    % if the linearised stability is to be computed, Graph should be the
    % list of edges.
    else
        if size(G,1) == size(G,2)
            [rows,cols,vals] = find(G);
            if sum(vals)==length(vals)
                Graph=[cols-1, rows-1 ones(size(cols))];
            else
                Graph=[cols-1, rows-1, vals];
            end
            clear rows cols vals;
        elseif size(G,2)==3
            Graph=G;
        else
            error('Wrong size for G: G should be a graph saved either as a list of edges (size(G)=[N,3] if weighted, size(G)=[N,2] if unweighted) or as an adjacency matrix (size(G)=[N,N])');
        end
        PARAMS.NbNodes = max(Graph(:,1))+1;
    end
else
    error('Please provide at least the graph to be partitioned. Type "help stability" for more information.');
end    

% Argument 2: T

if nargin > 1
    if (isvector(T) && isnumeric(T))
        Time=T;
    else
        error('The second argument should be a numerical vector. Type "help stability" for more information.');
    end
end

% Parallel computation: Initialize the number of cores if matlabpool is not
% yet running.
if PARAMS.ComputeParallel && (matlabpool('size') == 0)
    flag_matlabpool = true;
    matlabpool
end

%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%$                                          $%
%$      Computation of the stability        $%
%$                                          $%
%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%

% Display starting time
if verbose
    c = clock;
    disp(' ');
    disp(['   Partitioning of the graph started at ' datestr([2011 1 1 c(4) c(5) c(6)], 'HH:MM:SS') ' with the following parameters:']);
    disp(['      Stability function: ' func2str(StabilityFunction)]);
    disp(['      Compute the variation of information: ' int2str(PARAMS.ComputeVI)]);
    disp(['      Save the results in files: ' int2str(OutputFile)]);
    if OutputFile; disp(['      Prefix of the output files: ' prefix]); end
    disp(['      Number of Louvain iterations: ' int2str(PARAMS.NbLouvain)]);
    disp(['      Number of Louvain iterations used for the computation of VI: ' int2str(PARAMS.M)]);
    disp(['      Full stability: ' int2str(Full)]);
    disp(['      Check the input graph: ' int2str(Sanity)]);
    disp(['      Precision used: ' num2str(PARAMS.Precision)]);
    disp(['      Plot the results: ' int2str(plotStability)]);
    disp(['      Verbose mode: ' int2str(verbose)]);    
    disp(' ');
    tstart=tic;
end

if(Full)
    nr_edges = length(find(tril(Graph)~=0));
else
    nr_edges = length(Graph)/2;
end


% Initialisation
S = zeros(1, length(Time));
N = zeros(1, length(Time));
VI = zeros(1, length(Time));
C = zeros(PARAMS.NbNodes, length(Time));

if TextOutput
    mkdir(['Partitions_' prefix]);
end
if PARAMS.ComputeES
    ES = zeros(nr_edges, length(Time));
end
if OutputFile
    save(['Stability_' prefix '.mat'],'Time','S','N','VI','C');
end

if plotStability
    figure
end

if verbose
    step_prec=0;
end


% Loop over all Markov times
for t=1:length(Time)
    
    if verbose
        disp(['   Partitioning for Markov time = ' num2str(Time(t),'%10.6f') '...']);
    end
    
    
    
    [S(t), N(t), C(:,t), VI(t) VAROUT] = StabilityFunction(Graph, Time(t), PARAMS);
    if isfield(VAROUT,'precomputed')
        PARAMS.precomputed = VAROUT.precomputed;
        PARAMS.pi = VAROUT.pi;
        PARAMS.P = VAROUT.P;
    end
    
    if plotStability && t>1
        stability_plot(Time,t,S,N,VI,PARAMS.ComputeVI);
    end
    
    if TextOutput
        cd(['Partitions_' prefix]);
        dlmwrite(['Partition_' prefix '_' num2str(Time(t),'%10.6f') '.dat'],[[1:PARAMS.NbNodes]',C(:,t)],'delimiter','\t');
        cd ..;        
        dlmwrite(['Stability_' prefix '.stdout'],[Time(t), S(t), N(t), VI(t)],'-append', 'delimiter','\t')
    end   
    
    if OutputFile
        save(['Stability_' prefix '.mat'],'Time','S','N','VI','C','-append');
    end

    
    if verbose && 100*t/length(Time) >= step_prec+10        
        disp(' ');
        disp(['   Completed: ' num2str(round(100*t/length(Time)),10) '%']);
        remaining_time=toc(tstart)*(1-t/length(Time))/(t/length(Time));
        nb_hours = floor(remaining_time/3600);
        nb_min = floor((remaining_time - nb_hours*3600)/60);
        nb_sec = round(remaining_time - nb_hours*3600 - nb_min*60);
        disp(['   Estimated time remaining: ' datestr([2011  1 1 nb_hours nb_min nb_sec], 'HH:MM:SS')]);%num2str(nb_hours) ':' num2str(nb_min) ':' num2str(nb_sec)]);
        disp(' ');
        step_prec = floor(100*t/length(Time));
    end

end

if verbose
    c = clock;
    disp(' ');
    disp(['   Partitioning of the graph finished at ' datestr([2011 1 1 c(4) c(5) c(6)], 'HH:MM:SS')]);
    remaining_time=toc(tstart);
    nb_hours = floor(remaining_time/3600);
    nb_min = floor((remaining_time - nb_hours*3600)/60);
    nb_sec = round(remaining_time - nb_hours*3600 - nb_min*60);
    disp(['   Total time needed: ' datestr([2011 1 1 nb_hours nb_min nb_sec], 'HH:MM:SS')]);%num2str(nb_hours) ':' num2str(nb_min) ':' num2str(nb_sec)]);
end


if OutputFile
    save(['Stability_' prefix '.mat'],'Time','S','N','VI','C','-append');
end

if flag_matlabpool
    matlabpool close;
end



end

%------------------------------------------------------------------------------
function [StabilityFunction, OutputFile, prefix, Full, Sanity, plotStability, verbose, TextOutput,PARAMS] = parseinput(options,varargin)
% Parse the options from the command line

%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%
%$                                          $%
%$          Default parameters              $%
%$                                          $%
%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Global" options relevant for output and control flow
StabilityFunction = @louvain_FNL;    	        % Full stability with normalised laplacian is used by default
OutputFile = false;                             % No output file by default.
Laplacian = 'Normalised';                       % Default Laplacian
Full = true;                                    % If true, performs the full stability
Sanity = true;                                  % If true, performs the graph sanity checks
plotStability = false;                          % If true, plots the results of the stability, number of communities and variation of information vs Markov time.           
verbose = false;                                % Toggle verbose mode
prefix = '';                                    % Output prefix
TextOutput = false;                             % Toggles the text output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options stored in struct relevant for optimization etc.
PARAMS = struct;                                % create empty structure for storing parameters
PARAMS.precomputed = false;                     % Flag for precomputed transition matrix + stationary distribution
PARAMS.directed = false;                        % enables dealing with directed graphs
PARAMS.ComputeVI = true;                        % True if the variation of information should be computed
PARAMS.ComputeES = false;                       % True if edge statistics should be computed
PARAMS.ComputeParallel = false;                 % Toggles the computation in parallel
PARAMS.NbLouvain = 100;                         % Number of louvain optimisations at each Markov time
PARAMS.NbNodes = 0;                             % Total number of nodes;
PARAMS.Precision = 1e-9;                        % Threshold for stability and edges weigths
PARAMS.M = 100;                                 % Top M partitions among the L found by louvain are used to compute the variation of information
PARAMS.K = NaN;                                 % K stabilities value, only relevant for Ruelle random walk and k stabilities. K =-1 corresponds to the normalised Laplacian
PARAMS.teleport_tau = 0;                        % teleportation probability (only relevant for directed graphs)

% actual parsing begins here
attributes={'novi', 'l', 'm', 'out', 'linearised', 'nocheck', 'laplacian', 'prec', 'plot','v','t','p','k','directed','teleport','precomputed'};

if options > 0
    
    varargin = varargin{:}; % extract cell array input from varargin
    
    % test whether attribute-value pairs are specified, or fixed parameter order
    stringoptions = lower(varargin(cellfun('isclass',varargin,'char')));
    attributeindexesinoptionlist = ismember(stringoptions,attributes);
    newinputform = any(attributeindexesinoptionlist);
    if newinputform
        % parse values to functions parameters
        i = 1;
        while (i <= length(varargin))
            if strcmpi(varargin{i},'linearised')
                Full = false;
                i = i+1;
            elseif strcmpi(varargin{i},'directed')
                PARAMS.directed = true;
                i = i+1;
            elseif strcmpi(varargin{i},'novi')
                PARAMS.ComputeVI = false;
                i = i+1;
            elseif strcmpi(varargin{i},'nocheck')
                Sanity = false;
                i = i+1;
            elseif strcmpi(varargin{i},'plot')
                plotStability = true;
                i = i+1;
            elseif strcmpi(varargin{i},'v')
                verbose = true;
                i = i+1;
            elseif strcmpi(varargin{i},'precomputed')
                PARAMS.precomputed = true;
                i = i+1;
            elseif strcmpi(varargin{i},'p')
                if exist('matlabpool','file')
                    PARAMS.ComputeParallel = true;
                else
                    PARAMS.ComputeParallel = false;
                    warning('The Parallel Computing Toolbox of Matlab does not appear to be installed. Defaulting to single node computation...');
                end
                i = i+1;
            elseif strcmpi(varargin{i},'t')
                TextOutput = true;
                i = i+1;
            else
                %Check to make sure that there is a pair to go with
                %this argument.
                if length(varargin) < i + 1
                    error('MATLAB:stability:AttributeList', ...
                        'Attribute %s requires a matching value', varargin{i});
                elseif strcmpi(varargin{i},'laplacian')
                    if ischar(varargin{i+1})
                        Laplacian = varargin{i+1};
                    else
                        error('MATLAB:stability:laplacian',...
                            'Please provide a matching value for attribute laplacian. It must either be ''normalised'' or ''combinatorial''.');
                    end
                elseif strcmpi(varargin{i},'l')
                    if isnumeric(varargin{i+1})
                        PARAMS.NbLouvain = round(varargin{i+1});
                        PARAMS.M = round(varargin{i+1});
                    end
                elseif strcmpi(varargin{i},'prec')
                    if isnumeric(varargin{i+1})
                        PARAMS.Precision = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'m')
                    if isnumeric(varargin{i+1})
                        PARAMS.M = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'k')
                    if isnumeric(varargin{i+1})
                        PARAMS.K = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'teleport')
                    if isnumeric(varargin{i+1})
                        PARAMS.teleport_tau = varargin{i+1};
                    end
                elseif strcmpi(varargin{i},'out')
                    if ischar(varargin{i+1})
                        OutputFile = true;
                        prefix = varargin{i+1};
                    else
                        error('MATLAB:stability:out',...
                            'Please provide a matching value for attribute out. It must be a string.');
                    end
                else
                    error('MATLAB:stability:Attribute',...
                        'Invalid attribute tag: %s', varargin{i});
                end
                i = i+2;
            end
        end
    else 
        if ischar(varargin{1})
            error('MATLAB:stability:Attribute',...
                            'Invalid attribute tag: %s', varargin{1});
        else
            error('MATLAB:stability:Attribute',...
                            'Invalid attribute tag: %d', varargin{1});
        end
    end
end

TextOutput = TextOutput && OutputFile;

% Choose which type of stability is to be computed
if Full
    if strcmpi(Laplacian, 'combinatorial')
        StabilityFunction = @louvain_FCL;
    elseif strcmpi(Laplacian, 'normalised')
        StabilityFunction = @louvain_FNL;
    elseif strcmpi(Laplacian, 'corr_normalised')
        StabilityFunction = @louvain_FCNL;
    elseif strcmpi(Laplacian, 'mi_normalised')
        StabilityFunction = @louvain_MINL;
    elseif strcmpi(Laplacian, 'k_stability')
        StabilityFunction = @louvain_k_stability;
    elseif strcmpi(Laplacian, 'Ruelle_k_stability')
        StabilityFunction = @louvain_Ruelle_k_stability;
    else
        error('Please provide a valid matching value for attribute laplacian. It must either be ''normalised'' or ''combinatorial''.');
    end
else
    if strcmpi(Laplacian, 'combinatorial')
        StabilityFunction = @louvain_LCL;
    elseif strcmpi(Laplacian, 'normalised')
        StabilityFunction = @louvain_LNL;
    elseif strcmpi(Laplacian, 'corr_normalised')
        StabilityFunction = @louvain_CLNL;
    else
        error('Please provide a valid matching value for attribute laplacian.');
    end
end

end

%--------------------------------------------------------------------------
function shares = split_even(N,nr_threads)
%Function to compute an even split of N runs between t threads 
shares = ones(1,nr_threads)*floor(N/nr_threads);
shares(1) = shares(1)+ rem(N,nr_threads);
end

%TODO Implement edge statistics variant more fully, double check
%implementation of directed graph, parallel computations etc.
%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_FNL(Graph, time, PARAMS)
% Computes the full normalised stabilty

VAROUT =[]; % init varying outputs 

% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % directed case: M_ij >> from i to j
    if PARAMS.directed == true
        dout = sum(Graph,2);
        dangling = (dout==0);
        dout(dangling) = 1;
        Dout = sparse(diag(dout));
        clear dout;
        M = (1-PARAMS.teleport_tau)*Dout\Graph; % deterministic part of transition
        % teleportation according to arXiv:0812.1770
        M =	M + diag(PARAMS.teleport_tau + dangling.*(1-PARAMS.teleport_tau))...
            * ones(PARAMS.NbNodes)/PARAMS.NbNodes;
        
        clear Dout dangling
        [v lambda_all] = eigs(M'); % largest eigenvalue of transition matrix corresponds to stat.distribution.
        lambda = max(diag(lambda_all));
        v = v(:,diag(lambda_all) == lambda);
        v = abs(v);              % make sure eigenvector is positive
        clear lambda;
        % store results for future use
        VAROUT.precomputed = true;
        VAROUT.pi = v/sum(v);
        VAROUT.P = M;      
        % now compute exponential transition matrix
        solution = diag(v/sum(v))*expm(time* (M-eye(size(M))) );
        clear M v;
        % symmetrize solution
        solution = (solution+solution')/2;
        
        
        % undirected case
    else
        % Generate the matrix exponential
        trans=sparse(diag(    (sum(Graph)).^(-1)     ) * Graph);  %(stochastic) transition matrix        
        Lap=sparse(trans-eye(PARAMS.NbNodes));
        % store results for future use
        VAROUT.precomputed = true;
        VAROUT.P = trans;
        
        clear trans;
        exponential=sparse(expm(time.*Lap));
        clear Lap;
        
        PI=sparse((diag(sum(Graph)))/sum(sum(Graph)));  %diag matrix with stat distr   
        VAROUT.pi = diag(PI);   % store results for future use
        
        solution=sparse(PI*exponential);
        clear exponential;
        clear PI;
        
        
    end
    
    % stationary distribution and "transition matrix" have been computed before
else
    if PARAMS.directed == true
        solution = diag(PARAMS.pi)*expm(time* (PARAMS.P - eye(size(PARAMS.P))) );
        solution = (solution +solution')/2; % symetrization needed for directed case
    else
        solution = diag(PARAMS.pi)*expm(time* (PARAMS.P - eye(size(PARAMS.P))) );
    end
end



% prune out weights that are too small as defined by precision
solution=max(max(solution))*PARAMS.Precision*round(solution/(max(max(solution))*PARAMS.Precision));
[row,col,val] = find(solution);
clear solution
graph=[col-1,row-1,val];

% init some numbers
lnk = zeros(PARAMS.NbNodes, PARAMS.NbLouvain);
lnkS = zeros(PARAMS.NbLouvain,1);
nb_comm = zeros(PARAMS.NbLouvain,1);

% Optimize with Louvain NbLouvain times
if PARAMS.ComputeParallel
    nr_threads = matlabpool('size');
    stability = cell(1,nr_threads);
    nb_comm_temp = cell(1,nr_threads);
    communities = cell(1,nr_threads);
    shares = split_even(PARAMS.NbLouvain,nr_threads);
    precision = PARAMS.Precision;
    % computation in parallel with cell arrays
    parfor l=1:nr_threads;
        [stability{l}, nb_comm_temp{l}, communities{l}] = ...
            stability_louvain(graph, 1, shares(l), precision ,'normalised',randi(intmax));
    end    
    % assignements
    lnk = cat(2,communities{:});
    lnkS = cat(2,stability{:});
    nb_comm = cat(2,nb_comm_temp{:});
       
    % non parallel version
else
    [stability, nb_comm, communities] = ...
        stability_louvain(graph, 1, PARAMS.NbLouvain, PARAMS.Precision,'normalised',randi(intmax));
    lnk = communities;
    lnkS = stability;
end
clear communities stability nr_threads shares nb_comm_temp graph;


% Comment: maybe one should pick one of the best solutions at random,
% if two solutions have the same value;
index = find(lnkS==max(lnkS),1);

S = lnkS(index);
C = lnk(:,index);
N = nb_comm(index);



if PARAMS.ComputeVI && nnz(max(lnk)==PARAMS.NbNodes-1)~=PARAMS.NbLouvain && nnz(max(lnk)==0)~=PARAMS.NbLouvain
    VI = computeRobustness(lnk, lnkS, PARAMS.M,PARAMS.ComputeParallel);
else
    VI=0; 
end

clear lnk;

end

%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_FCL(Graph, time, PARAMS)
% Computes the full combinatorial stability

VAROUT =[]; % init varying outputs 

% Directed part so far "standard implementation" with out degree Laplacian,
% check for bugs etc...

% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % directed case: M_ij >> from i to j
    if PARAMS.directed == true
        dout = sum(Graph,2);
        dangling = (dout==0);
        av_degree = sum(dout)/PARAMS.NbNodes; 
        clear dout;
        M = (1-PARAMS.teleport_tau)*Graph; % deterministic part of transition
        % teleportation according to arXiv:0812.1770
        M =	M + av_degree*diag(PARAMS.teleport_tau + dangling.*(1-PARAMS.teleport_tau))...
            * ones(PARAMS.NbNodes)/PARAMS.NbNodes;
        % dout of new Graph
        dout = sum(M,2);
        Dout = sparse(diag(dout));
        av_degree = sum(dout)/PARAMS.NbNodes;         
        % out degree Laplacian normalized by new average degree
        Lap = (Dout -Graph)/av_degree;
        clear Dout dangling
        [v lambda_all] = eigs(Lap',5,'SM'); % zero eigenvalue of Laplacian matrix corresponds to stat.distribution.
        lambda = min(diag(lambda_all));
        v = v(:,diag(lambda_all) == lambda);
        v = abs(v);              % make sure eigenvector is positive
        clear lambda;
        % store results for future use
        VAROUT.precomputed = true;
        VAROUT.pi = v/sum(v);
        VAROUT.P = Lap;      
        % now compute exponential transition matrix
        solution = diag(v/sum(v))*expm(-time*Lap );
        clear M v Lap;
        % symmetrize solution
        solution = (solution+solution')/2;
        
        
        % undirected case
    else
        % standard Laplacian and average degree
        av_degree = sum(sum(Graph))/PARAMS.NbNodes;
        Lap=  sparse(Graph-diag(sum(Graph)));
        % store results for future use
        VAROUT.precomputed = true;
        VAROUT.P = Lap/av_degree;
        
 
        % Generate the matrix exponential
        exponential=sparse(expm(-time.*Lap/av_degree));
        clear Lap;
        
        PI=sparse(eye(PARAMS.NbNodes)/PARAMS.NbNodes);  %diag matrix with stat distr   
        VAROUT.pi = diag(PI);   % store results for future use
        
        solution=sparse(PI*exponential);
        clear exponential;
        clear PI;
        
        
    end
    
    % stationary distribution and "transition matrix" have been computed before
else
    if PARAMS.directed == true
        solution = diag(PARAMS.pi)*expm(-time* PARAMS.P);
        solution = (solution +solution')/2; % symetrization needed for directed case
    else
        solution = diag(PARAMS.pi)*expm(-time* PARAMS.P);
    end
end



% prune out weights that are too small as defined by precision
solution=max(max(solution))*PARAMS.Precision*round(solution/(max(max(solution))*PARAMS.Precision));
[row,col,val] = find(solution);
clear solution
graph=[col-1,row-1,val];

% init some numbers
lnk = zeros(PARAMS.NbNodes, PARAMS.NbLouvain);
lnkS = zeros(PARAMS.NbLouvain,1);
nb_comm = zeros(PARAMS.NbLouvain,1);

% Optimize with Louvain NbLouvain times
if PARAMS.ComputeParallel
    nr_threads = matlabpool('size');
    stability = cell(1,nr_threads);
    nb_comm_temp = cell(1,nr_threads);
    communities = cell(1,nr_threads);
    shares = split_even(PARAMS.NbLouvain,nr_threads);
    precision = PARAMS.Precision;
    % computation in parallel with cell arrays
    parfor l=1:nr_threads;
        [stability{l}, nb_comm_temp{l}, communities{l}] = ...
            stability_louvain(graph, 1, shares(l), precision ,'normalised',randi(intmax));
    end    
    % assignements
    lnk = cat(2,communities{:});
    lnkS = cat(2,stability{:});
    nb_comm = cat(2,nb_comm_temp{:});
       
    % non parallel version
else
    [stability, nb_comm, communities] = ...
        stability_louvain(graph, 1, PARAMS.NbLouvain, PARAMS.Precision,'normalised',randi(intmax));
    lnk = communities;
    lnkS = stability;
end
clear communities stability nr_threads shares nb_comm_temp graph;


% Comment: maybe one should pick one of the best solutions at random,
% if two solutions have the same value;
index = find(lnkS==max(lnkS),1);

S = lnkS(index);
C = lnk(:,index);
N = nb_comm(index);



if PARAMS.ComputeVI && nnz(max(lnk)==PARAMS.NbNodes-1)~=PARAMS.NbLouvain && nnz(max(lnk)==0)~=PARAMS.NbLouvain
    VI = computeRobustness(lnk, lnkS, PARAMS.M,PARAMS.ComputeParallel);
else
    VI=0; 
end

clear lnk;

end
%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_FCNL(Graph, time, PARAMS)
% Computes the full normalised corr stabilty

VAROUT =[]; % init varying outputs 
%TODO adjust below code properly to work with paremters struct, so far just
% copy
ComputeES = PARAMS.ComputeES;
ComputeVI = PARAMS.ComputeVI  ;
precision = PARAMS.Precision;
NbLouvain = PARAMS.NbLouvain;
M = PARAMS.M ;
NbNodes = PARAMS.NbNodes;
ComputeParallel = PARAMS.ComputeParallel;

% Generate the matrix exponential
diagdeg=sparse((diag(sum(Graph)))/sum(sum(Graph)));  %diag matrix with stat distr
trans=sparse(diag(    (sum(Graph)).^(-1)     ) * Graph);  %(stochastic) transition matrix
clear Graph;
Lap=sparse(trans-eye(NbNodes));
clear trans;
exponential=sparse(expm(time.*Lap));
clear Lap;
solution=sparse(diagdeg*exponential);
clear exponential;
clear diagdeg;
solution=max(max(solution))*precision*round(solution/(max(max(solution))*precision));
clear exponential;
clear diagdeg;
[row,col,val] = find(solution);
clear solution
graph=[col-1,row-1,val];

% Optimize louvain NbLouvain times
[stability, nb_comm, communities] = stability_louvain(graph, 1, NbLouvain, precision,'corr_normalised');
lnk = communities;
lnkS = stability;
% Comment: maybe one should pick one of the best solutions at random,
% if two solutions have the same value;
index = find(stability==max(stability),1);

S = stability(index);
C = communities(:,index);
N = nb_comm(index);

clear communities;
clear graph;

if ComputeVI && nnz(max(lnk)==NbNodes-1)~=NbLouvain && nnz(max(lnk)==0)~=NbLouvain
     VI = computeRobustness(lnk, lnkS, M, ComputeParallel);
else
    VI=0;
end

clear lnk;

end

%THIS FUNCTION IS A TESTBED ATM!!!!!!!!!!!! BE AWARE!
%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_MINL(Graph, time, PARAMS)
% Computes the full normalised mutual information stabilty

VAROUT =[]; % init varying outputs 
%TODO adjust below code properly to work with paremters struct, so far just
% copy
ComputeES = PARAMS.ComputeES;
ComputeVI = PARAMS.ComputeVI  ;
precision = PARAMS.Precision;
NbLouvain = PARAMS.NbLouvain;
M = PARAMS.M ;
NbNodes = PARAMS.NbNodes;
ComputeParallel = PARAMS.ComputeParallel;

%TODODODODODODOO
% Generate the matrix exponential
diagdeg=sparse((diag(sum(Graph)))/sum(sum(Graph)));  %diag matrix with stat distr
trans=sparse(diag(    (sum(Graph)).^(-1)     ) * Graph);  %(stochastic) transition matrix

Lap=sparse(trans-eye(NbNodes));
clear trans;
exponential=sparse(expm(time.*Lap));
%solution=sparse(Lap*diagdeg*exponential);% changes in here
solution=sparse(diagdeg*exponential);% changes in here
clear Lap;
M_null = ones(size(sum(Graph)))'*sum(Graph)/sum(Graph(:));
pi = sum(Graph)/sum(sum(Graph));
solution = (solution+pi'*pi).*log2(solution./(pi'*pi));
% solution = solution.*log2(solution./(pi'*pi));% +(eye(size(M_null))-M_null)*diagdeg*expm(-(eye(size(M_null))-M_null)*time);
%solution = solution - diagdeg*expm(-(eye(size(M_null))-M_null)*time*.9);
solution = createNormalizedLyapunovStabilityMatrix(Graph,time);
solution = (solution+solution')/2;
clear exponential;
clear diagdeg;
%solution=max(max(solution))*precision*round(solution/(max(max(solution))*precision));
clear exponential;
clear diagdeg;
% null_model = sum(Graph)'*sum(Graph)/(sum(Graph(:))^2);
% solution = solution.*log2(solution./(null_model));
% solution(isnan(solution))=0;
clear Graph null_model;
[row,col,val] = find(solution);
clear solution
graph=[col-1,row-1,val];

% Optimize louvain NbLouvain times
[stability, nb_comm, communities] = stability_louvain(graph, time, NbLouvain, precision,'mutual_information');
lnk = communities;
lnkS = stability;
% Comment: maybe one should pick one of the best solutions at random,
% if two solutions have the same value;
index = find(stability==max(stability),1);

S = stability(index);
C = communities(:,index);
N = nb_comm(index);

clear communities;
clear graph;

if ComputeVI && nnz(max(lnk)==NbNodes-1)~=NbLouvain && nnz(max(lnk)==0)~=NbLouvain
     VI = computeRobustness(lnk, lnkS, M,ComputeParallel);
else
    VI=0;
end

clear lnk;

end


%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_k_stability(Graph, time, PARAMS)
%Computes the full k-Laplacian stabilty

VAROUT =[]; % init varying outputs 
if(isnan(PARAMS.K))
    error('Please provide a value for K if using k-Laplacian stability');
elseif(PARAMS.K == -1)
    [S, N, C, VI] = louvain_FCL(Graph, time, PARAMS);
elseif(PARAMS.K == 0)
    [S, N, C, VI] = louvain_FNL(Graph, time, PARAMS);
else
    k = PARAMS.K;
    precision = PARAMS.Precision;
    
    % Generate the matrix exponential
    D=sparse(diag(sum(Graph)));  %diag matrix with stat distr
    M=sparse(diag(  (sum(Graph)).^(-1) ) * Graph);  %(stochastic) transition matrix
    d_mk_av = mean( diag(D).^(-k) ); % vector with mean of degree to the minus k-th power 
    Lap=diag(diag(D).^(-k))  \ sparse(M-eye(PARAMS.NbNodes)); % k-Laplacian, actually minus L_k
    clear M;
    exponential=sparse(expm(time*Lap/d_mk_av));
    clear Lap;
    clear d_mk_av;
    solution=sparse(diag(diag(D).^(k+1))*exponential);
    clear exponential;
    clear D;
    clear k;
    solution=max(max(solution))*precision*round(solution/(max(max(solution))*precision));
    clear exponential;
    [row,col,val] = find(solution);
    clear solution
    graph=[col-1,row-1,val];
    
    % Optimize louvain NbLouvain times
    [stability, nb_comm, communities] = stability_louvain(graph, 1, PARAMS.NbLouvain, precision,'normalised');
    lnk = communities;
    lnkS = stability;
    % Comment: maybe one should pick one of the best solutions at random,
    % if two solutions have the same value;
    index = find(stability==max(stability),1);
    
    S = stability(index);
    C = communities(:,index);
    N = nb_comm(index);
    
    clear communities;
    clear graph;
    
    if PARAMS.ComputeVI% && nnz(max(lnk)==NbNodes-1)~=NbLouvain && nnz(max(lnk)==0)~=NbLouvain
        %[VI, nr_cores, cores, edge_statistics] = findCoreAndPeriphery(Graph,lnk);
        VI = computeRobustness(lnk, lnkS, PARAMS.M,PARAMS.ComputeParallel);
    else
        VI=0;
    end
    
    clear lnk;
end
end

%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_Ruelle_k_stability(Graph, time, PARAMS)
% Computes the full k-Laplacian Ruelle stabilty
VAROUT =[]; % init varying outputs 
if(isnan(PARAMS.K))
    error('Please provide a value for K if using k-Ruelle stability');
else
    k = PARAMS.K;
    precision = PARAMS.Precision;
    
    % largest eigenvalue of graph adjacency and corresponding EV
    [v lambda_all] = eigs(Graph);         % be careful with eigs as results are ordered according to magnitude!
    lambda = max(diag(lambda_all));
    v = v(:,diag(lambda_all) == lambda);
    LAMBDA_v =diag(abs(v));           
    clear v;
    M_R = 1/lambda * (LAMBDA_v\Graph*LAMBDA_v); % Ruelle transition matrix;
    clear lambda lambda_all;
    
    % Generate the matrix exponential
    v_mk_av = mean( diag(LAMBDA_v).^(-k) ); % mean of v to the minus k-th power 
    Lap=diag(diag(LAMBDA_v).^(-k)) \ sparse(M_R-eye(PARAMS.NbNodes)); % Ruelle-k-Laplacian, actually negative of it
    clear M_R;
    exponential=sparse(expm(time*Lap/v_mk_av));
    clear Lap;
    clear v_mk_av;
    pi = diag(LAMBDA_v).^(k+2)/sum(diag(LAMBDA_v).^(k+2));
    PI= sparse(diag(pi));
    clear pi LAMBDA_v;
    solution=PI*exponential;    
    clear PI exponential k;
    solution=max(max(solution))*precision*round(solution/(max(max(solution))*precision));
    clear exponential;
    [row,col,val] = find(solution);
    clear solution
    
    % adjust range of values, important as otherwise val tend to become too small..
%      mval = mean(val);
%      val = val/mval; clear mval
    
    graph=[col-1,row-1,val];
        
    % Optimize louvain NbLouvain times
    [stability, nb_comm, communities] = stability_louvain(graph, 1, PARAMS.NbLouvain, precision,'normalised');
    lnk = communities;
    lnkS = stability;
    % Comment: maybe one should pick one of the best solutions at random,
    % if two solutions have the same value;
    index = find(stability==max(stability),1);
    
    S = stability(index);
    C = communities(:,index);
    N = nb_comm(index);
    
    clear communities;
    clear graph;
    
    if PARAMS.ComputeVI% && nnz(max(lnk)==NbNodes-1)~=NbLouvain && nnz(max(lnk)==0)~=NbLouvain
        %[VI, nr_cores, cores, edge_statistics] = findCoreAndPeriphery(Graph,lnk);
        VI = computeRobustness(lnk, lnkS, PARAMS.M,PARAMS.ComputeParallel);
    else
        VI=0;
    end
    
    clear lnk;
end
end

%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_LCL(Graph, time, PARAMS)
VAROUT =[]; % init varying outputs 
%TODO adjust below code properly to work with paremters struct, so far just
% copy
ComputeES = PARAMS.ComputeES;
ComputeVI = PARAMS.ComputeVI  ;
precision = PARAMS.Precision;
NbLouvain = PARAMS.NbLouvain;
M = PARAMS.M ;
NbNodes = PARAMS.NbNodes;
ComputeParallel = PARAMS.ComputeParallel;

% Optimize louvain NbLouvain times
[stability, nb_comm, communities] = stability_louvain(Graph, time, NbLouvain, precision,'combinatorial');
lnk = communities;
lnkS = stability;
% Comment: maybe one should pick one of the best solutions at random,
% if two solutions have the same value;
index = find(stability==max(stability),1);

S = stability(index);
C = communities(:,index);
N = nb_comm(index);


clear communities;
clear Graph;

if ComputeVI && nnz(max(lnk)==NbNodes-1)~=NbLouvain && nnz(max(lnk)==0)~=NbLouvain
     VI = computeRobustness(lnk, lnkS, M,ComputeParallel);
else
    VI = 0;
end

clear lnk;

end
%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_LNL(Graph, time, PARAMS)
VAROUT =[]; % init varying outputs 
%TODO adjust below code properly to work with paremters struct, so far just
% copy
ComputeES = PARAMS.ComputeES;
ComputeVI = PARAMS.ComputeVI  ;
precision = PARAMS.Precision;
NbLouvain = PARAMS.NbLouvain;
M = PARAMS.M ;
NbNodes = PARAMS.NbNodes;
ComputeParallel = PARAMS.ComputeParallel;

% Optimize louvain NbLouvain times
[stability, nb_comm, communities] = stability_louvain(Graph, time, NbLouvain, precision,'normalised');
lnk = communities;
lnkS = stability;
% Comment: maybe one should pick one of the best solutions at random,
% if two solutions have the same value;
index = find(stability==max(stability),1);

S = stability(index);
C = communities(:,index);
N = nb_comm(index);



clear communities;
clear Graph;

if ComputeVI && nnz(max(lnk)==NbNodes-1)~=NbLouvain && nnz(max(lnk)==0)~=NbLouvain
    VI = computeRobustness(lnk,lnkS, M,ComputeParallel);
else
    VI = 0;
end

clear lnk;

end
%------------------------------------------------------------------------------
function [S, N, C, VI, VAROUT] = louvain_CLNL(Graph, PARAMS)
VAROUT =[]; % init varying outputs 
%TODO adjust below code properly to work with paremters struct, so far just
% copy
ComputeES = PARAMS.ComputeES;
ComputeVI = PARAMS.ComputeVI  ;
precision = PARAMS.Precision;
NbLouvain = PARAMS.NbLouvain;
M = PARAMS.M ;
NbNodes = PARAMS.NbNodes;
ComputeParallel = PARAMS.ComputeParallel;

% Optimize louvain NbLouvain times
[stability, nb_comm, communities] = stability_louvain(Graph, time, NbLouvain, precision,'corr_normalised');
lnk = communities;
lnkS = stability;
% Comment: maybe one should pick one of the best solutions at random,
% if two solutions have the same value;
index = find(stability==max(stability),1);

S = stability(index);
C = communities(:,index);
N = nb_comm(index);



clear communities;
clear Graph;

if ComputeVI && nnz(max(lnk)==NbNodes-1)~=NbLouvain && nnz(max(lnk)==0)~=NbLouvain
    VI = computeRobustness(lnk,lnkS, M,ComputeParallel);
else
    VI = 0;
end

clear lnk;

end
%------------------------------------------------------------------------------
function Graph = check(Graph, verbose, PARAMS)
% Check that the graph is properly encoded.
    if verbose
        disp(' ');
        disp('   Graph sanity check...');
    end
    
    % Initialisation of Graph properties
    edgelist = false;
    unweighted = false;
    
    if size(Graph,2) < 2 || size(Graph,1) < 3
        error(['The size of the graph is [' num2str(size(Graph,1)) ',' num2str(size(Graph,2)) ']. Please check that it has been correclty encoded.'])
    end       
    
    if size(Graph,2) ~= size(Graph,1)
        if size(Graph,2) ~=2 && size(Graph,2) ~=3
            error('Wrong size for G: G should be a graph saved either as a list of edges (size(G)=[N,3] if weighted, size(G)=[N,2] if unweighted) or as an adjacency matrix (size(G)=[N,N])');
        end
        edgelist = true;
        if size(Graph,2) == 2
            unweighted = true;
        end
    end
    
    % Check nodes numbering and convert edgelist into adjacency matrix
    if edgelist
        if min(min(Graph(:,1:2))) ~=0
            warning('The numbering of the nodes in the graph should always start with zero.');
            old_node_1 = min(min(Graph(:,1:2)));
            Graph(:,1)=Graph(:,1)-old_node_1;
            Graph(:,2)=Graph(:,2)-old_node_1;
        end
        if unweighted == false
            Graph=sparse(Graph(:,1)+1,Graph(:,2)+1,Graph(:,3));
        else
            Graph=sparse(Graph(:,1)+1,Graph(:,2)+1,ones(length(Graph(:,1)),1));
        end
    end

    % Check that graph contains just numbers
    if any(any(~isnumeric(Graph)))
	error('The graph provided contains elements which are not numbers (isnumeric == false). Please check your graph, and try again.');
    end
        
    % Check symmetry of the adjacency matrix if graph is not directed
    if PARAMS.directed == false
    	if size(Graph,1) ~= size(Graph,2)
        	error('The graph provided is a directed graph. Specify the correct options or correct your graph');
    	end
	if any(any(Graph~=Graph'))
		if nnz(triu(Graph,1))>0 && nnz(tril(Graph,-1))>0
		    error('The graph provided is a directed graph.');
		else
		    warning('Adjacency matrix A of provided graph is triangular -- symmetrizing A = A + A^T');
		    Graph=Graph+Graph';
		end
	end
    end

    % Check for isolated nodes
    if nnz(sum(Graph))~=size(Graph,2)
        warning('There are isolated nodes in the graph');
    end
    
    % Check for disconnected components
    if exist('graphconncomp','file') == 2
        nbcomp=graphconncomp(sparse(Graph));
        if nbcomp>1
            warning(['There are ' num2str(nbcomp) ' not strongly connected components in the graph. If your graph is directed please be aware of the teleportation settings.']);
        end
    end
    
    % Return Graph to its original form
    if edgelist
        [row, col, val] = find(Graph);
        if unweighted
            Graph=[col-1, row-1];
        else
            Graph=[col-1, row-1, val];
        end
    else
        Graph=sparse(Graph);
    end
end
%------------------------------------------------------------------------------
function VI = computeRobustness(lnk, lnkS, M,ComputeParallel)

% Parameter

[~,i] = sort(lnkS);
lnk=lnk(:,i);
lnk=lnk(:,end-M+1:end);
VI = varinfo(lnk',ComputeParallel);
clear i;
end

%------------------------------------------------------------------------------
function [] = stability_plot(Time,t,S,N,VI,ComputeVI)

if ComputeVI
    subplot(2,1,1), ax=plotyy(Time(1:t),N(1:t),Time(N>1),S(N>1));
else
    ax=plotyy(Time(1:t),N(1:t),Time(N>1),S(N>1));
end
xlabel('Markov time');
set(ax(1),'YScale','log');%,'YTick',10.^[0:1:ceil(log10(max(N)))], 'YMinorTick', 'on');
set(ax(2),'YScale','log');%,'YTick',10.^[floor(log10(min(S(S>0)))):1:ceil(log10(max(S(S>0)))) ], 'YMinorTick', 'on');
set(ax(1),'YTickMode','auto','YTickLabelMode','auto','YMinorGrid','on');
set(ax(2),'YTickMode','auto','YTickLabelMode','auto','YMinorGrid','on');
set(get(ax(1),'Ylabel'),'String','Number of communities');
set(get(ax(2),'Ylabel'),'String','Stability');
% TODO compute limit, as corr stability has different normalization this
% has to be recalculated..
limits = [10^floor(log10(min(S(N>1)))) 1];
set(ax(1),'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YLim', [1 10^ceil(log10(max(N)))], 'XScale','log','XMinorGrid','on');
set(ax(2),'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YLim', limits, 'XScale','log');
ylabel('Number of communities');
if ComputeVI 
    subplot(2,1,2), semilogx(Time(1:t),VI(1:t));
    set(gca, 'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YMinorGrid','on','XMinorGrid','on');
    if max(VI)>0
        set(gca,'YLim', [0 max(VI)*1.1]);
    end
    xlabel('Markov time');
    ylabel('Variation of information');
end
drawnow;
end
