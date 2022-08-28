function [BestF,BestI,RunData,Cache] = mads(Problem,iterate0,varargin)
%MADS  Solver for nonlinear and mixed variable constrained optimization
%
%   Syntax:
%      [BESTF,BESTI,RUNDATA]       = mads(PROBLEM,ITERATE0)
%      [BESTF,BESTI,RUNDATA]       = mads(PROBLEM,ITERATE0,OPTIONS)
%      [BESTF,BESTI,RUNDATA]       = mads(PROBLEM,ITERATE0,OPTIONS,RUNDATA)
%      [BESTF,BESTI,RUNDATA,CACHE] = mads(PROBLEM,ITERATE0)
%      [BESTF,BESTI,RUNDATA,CACHE] = mads(PROBLEM,ITERATE0,OPTIONS)
%      [BESTF,BESTI,RUNDATA,CACHE] = mads(PROBLEM,ITERATE0,OPTIONS,RUNDATA)
%
%   Description:
%      MADS is an optimization algorithm for constrained nonlinear programming
%      (NLP) and mixed variable programming (MVP) problems.  It employs the
%      derivative-free class of mesh-adaptive direct search (MADS) filter 
%      algorithms, which are a generalization of the class of generalized 
%      pattern search (GPS) methods.  Derivatives are not necessary to run the
%      algorithm, but if available can be used to make the algorithm faster and
%      more efficient in finding a solution.
%
%      All the input and output variables are structures, each of which
%      contains many parameters.  The variables BESTF and BESTI, contain the
%      data describing the best feasible point and least infeasible point found
%      by the algorithm.  RUNDATA contains statistics on the MADS run, such as
%      final mesh size, number of iterations, number of function and gradient
%      evaluations, and CPU time used.  CACHE contains data associated with
%      every point evaluated during the run.
%
%      The input variable PROBLEM contains variables that describe the
%      optimization problem to be solved.  In particular, it contains the names
%      of all the files associated with the problem, such as the functions
%      file, the linear constraints file, the discrete neighbors file for MVP
%      problems, and the initial points file.  ITERATE0 is a structure
%      describing all the initial iterates.  OPTIONS contains all user options
%      and settings to be used during the run.  If not used as an input
%      argument, MADS uses default settings, as specified in the MADS_DEFAULTS
%      function.  Using RUNDATA as an input variable is only done within the
%      NOMADm GUI setup, when resuming a previously stopped run.  It is highly
%      discouraged when running in batch mode.
%
%      The user is referred to the User's Guide for more detailed help in
%      understanding the input and output variables.
%
%   See also MADS_BATCH, MADS_DEFAULTS, NOMADM

%*******************************************************************************
%   Copyright (c) 2001-2016 by Mark A. Abramson
%
%   This file is part of the NOMADm software package.
%
%   NOMADm is free software; you can redistribute it and/or modify it under the
%   terms of the GNU General Public License as published by the Free Software
%   Foundation; either version 3 of the License, or (at your option) any later
%   version.
%
%   NOMADm is distributed in the hope that it will be useful, but WITHOUT ANY
%   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
%   details.  A copy of the GNU General Public License is available at
%   http://www.gnu.org/licenses/gpl-3.0.en.html.
%
% ------------------------------------------------------------------------------
%   Originally created, 2001.
%   Last modified, 11 January 2016
%
%   Author information:
%   Mark A. Abramson, LtCol (ret.), USAF, PhD
%   Abramson.Mark@gmail.com
%*******************************************************************************

%*******************************************************************************
% mads: Runs the Mesh-Adpative Direct Search (MADS) Progressive Barrier
%       Algorithm for solving constrained Mixed Variable Programming (MVP)
%       problems.
% ------------------------------------------------------------------------------
% CALLED FUNCTIONS:
%  processInput           = Initialize and validate input data and CACHE
%    mads_defaults        =   Retrieves MADS parameter default settings
%    checkErrors          =   Error check the input data
%    createCache          =   Create and initialize the Cache
%    createFilter         =   Create and initialize a filter
%    makeFeasible         =   Convert infeasible starting point to feasible
%  processOutput          = Process output and delete temporary files/variables
%    plotHistory          =   Plot objective value vs number of function evals
%  closeWorkspace         = Close appdata and delete temporary files
%  search                 = Perform MADS Search step
%    lhs                  =   Perform Latin Hypercube Search on the MADS mesh
%    es                   =   Perform Evolutionary Strategy
%      setPenaltyTerm     =     Set penalty parameter for constrained problem
%    getDataSites         =   Generate data sites used for forming surrogates
%    daceSurrogate        =   Construct a DACE surrogate
%    customSurrogate      =   Construct a custom surrogate
%      recalSurrogate     =     Construct or recalibrate a surrogate
%    optimizeSurrogate    =   Solve a surrogate optimization problem
%  poll                   = Perform MADS Poll step
%    standDirections      =   Retrieve standard mesh directions
%      activeConstraints  =     Identify e-active (linear) constraints
%      getTangentCone     =     Compute generators for the tangent cone
%        removeRedundancy =       Remove redundant active constraints
%      getAlternateCone   =     Compute generators of an alternative cone
%      scaleDirections    =     Scale the directions used in the Poll step
%    gradDirections       =   Retrieve gradient-pruned mesh directions
%      dVector            =     Compute the appropriate descent direction
%    madsDirections       =   Retrieve MADS mesh directions
%    getPollOrder         =   Set the order in which iterates are polled
%  mvpPoll                = Perform MADS discrete neighbor & Extended Poll step
%    < User N >           =   User-defined problem discrete neighbor file
%  update                 = Update MADS parameters and run statistics
%    getPollCenter        =   Determine the next poll center
%  updateOmega            = Update Omega and scale factors
%    getScaleFactors      =   Compute scale factors, based on initial data
%    < User O >           =   User-defined problem Omega file
%  evalPointSet           = Evaluate functions at multiple trial points
%    inOmega              =   Test if feasible w.r.t. linear constraints
%    isCacheHit           =   Test if point is already in the Cache
%    evalFunction         =   Evaluate problem functions at a given point
%      < User F >         =     User-defined problem Function file
%    updateCache          =   Update the Cache with the current iterate
%    updateFilter         =   Update the filter and solutions vectors
%      dominates          =     Tests if one iterate dominates another
%      plotFilter         =     Update the current real-time filter plot
%  terminate              = Tests if termination criteria are satisfied
%  runLP                  = Run a linear programming problem
%  cat2cont               = Convert categorical variables to continuous
%  dealIterate            = Deal matrix columns into vector of iterates
%  myeval                 = Evaluate variable of unknown type
%  default_N              = Default discrete neighbor file for MVP problems
%  cmaes_mod              = Evolutionary Strategy solver
%    xintobounds          =   Move infeasible point to feasible region boundary
%    myprctile            =   Return data value of specified percentile
% ------------------------------------------------------------------------------
%  linprog (Optim_TB)     = Optimization Toolbox LP solver
%  fmincon (Optim_TB)     = Optimization Toolbox constrained NLP solver
%  fminunc (Optim_TB)     = Optimization Toolbox unconstrained NLP solver
% ------------------------------------------------------------------------------
% USER FUNCTION FILES needed for Optimization Problem: example
%   example.m             = returns f(x,p), c(x,p), f'(x,p), c'(x,p)
%   example_Omega.m       = returns A,l,u,plist, given parameter p
%   example_X.m           = returns yes/no feasibility of any closed constraints
%   example_x0.m          = returns initial iterate(s)
%   example_N.m           = returns the set N of neighbor points 
%   example_Param.m       = returns a structure Param of user parameters
%   <Search file>         = optional function specifying Search procedure
% ------------------------------------------------------------------------------
% VARIABLES (only for the mads function):
%  BestF                  = final best feasible solution found
%  BestI                  = final least infeasible solution found
%  RunData                = statistics measuring algorithm performance
%    .delta               =   poll size
%    .nIter               =   cumulative number of iterations
%    .nFunc               =   cumulative number of function evaluations
%    .time                =   cumulative amount of CPU time expired
%    .nFails              =   number of consecutive Poll failures
%    .pcenter             =   current poll center
%    .stopRun             =   flag for stopping run at the next iteration
%  Cache                  = structure containing previously evaluated iterates
%    .Filter              =   structure containing filter data
%  Problem                = structure containing optimization problem data
%    .nameCache           =   name of the base workspace Cache variable
%    .typeProblem         =   string indicating "truth" or "surrogate"
%    .isMVP               =   flag indicating a mixed variable problem
%  iterate0               = initial iterate
%  Options                = structure containing SEARCH and POLL parameters
%    .Term                =   criteria for terminating the MADS function
%  success                = current iteration success or failure
%*******************************************************************************

% Initialize and verify validity of program variables
try
   [Problem,Options,RunData] = processInput(Problem,iterate0,varargin{:});

   % Main loop: Search, Poll, and Update parameters
   fmt1 =  '%sBest f-Value =  %25.16f %25.16f %25.16f\n';
   fmt2 = ['%sBest x-Value = [',repmat('%25.16f ',1,RunData.pcenter.n),']\n'];
   while ~terminate(Options.Term,RunData)
      if Options.debug   
         pLabel = sprintf('%s (%3d):  ',Problem.typeProblem,RunData.nIter+1);
         fprintf(fmt1,pLabel,RunData.pcenter.f,RunData.pcenter.h,RunData.delta);
         fprintf(fmt2,pLabel,RunData.pcenter.x');
      end

      if Options.debug >= 1, disp(' '); end
      if Options.debug >= 2
         disp([pLabel,'Mesh size = ',num2str(RunData.delta,Options.digits)]);
      end

      % Begin SEARCH step
      Cache = getappdata(0,Problem.nameCache);
      if Options.debug >= 1
         fprintf('%sBegin Search\n',pLabel);
      end
      [success,TempFilter,termflag] = search(Problem,Options,RunData, ...
                                             Cache.Filter);
      if termflag, RunData.stopRun = 1; continue, end
      if Options.debug >= 1
         fprintf('%sEnd   Search (success = %d)\n',pLabel,success);
      end
      Cache = getappdata(0,Problem.nameCache);
      Cache.Filter = TempFilter;
      setappdata(0,Problem.nameCache,Cache);

      % Begin POLL step
      if success < 0
         RunData.nPoll = RunData.nPoll + 1;
         if Options.debug >= 1
            fprintf('%sBegin Poll\n',pLabel);
         end
         [success,TempFilter,RunData] = poll('P',Problem,Options,RunData,...
                                             RunData.pcenter,Cache.Filter);
         if Options.debug >= 1
            fprintf('%sEnd   Poll (success = %d)\n',pLabel,success);
         end
         Cache = getappdata(0,Problem.nameCache);
         Cache.Filter = TempFilter;
         setappdata(0,Problem.nameCache,Cache);
      end

      % Begin EXTENDED POLL step (if MVP problem)
      if success < 0 && Problem.isMVP
         if Options.debug >= 1
            fprintf('%sBegin MV Poll\n',pLabel);
         end
         success = mvpPoll(Problem,Options,RunData);
         if Options.debug >= 1
            fprintf('%sEnd   MV Poll (success = %d)\n',pLabel,success);
         end
      end

      % Update MADS parameters
      [Problem,RunData] = update(Problem,Options,RunData,success);
   end
   if Options.debug   
      pLabel = sprintf('%s (%3d):  ',Problem.typeProblem,RunData.nIter);
      fprintf(fmt1,pLabel,RunData.pcenter.f,RunData.pcenter.h,RunData.delta);
      fprintf(fmt2,pLabel,RunData.pcenter.x');
   end

   % Post-processing
   [BestF,BestI] = processOutput(Problem,Options,RunData);
   Cache = closeWorkspace(Problem);

catch exception
   Cache = closeWorkspace(Problem);
   rethrow2(exception);
end
return


%*******************************************************************************
% FUNCTIONS FOR PROCESSING INPUT DATA
%*******************************************************************************
%*******************************************************************************
% processInput:  Initializes, validates, and processes input data
% ------------------------------------------------------------------------------
% Called by: mads
% Calls:     madsDefaults, checkErrors,  updateOmega,   createCache,
%            evalPointSet, makeFeasible, getPollCenter, update
% VARIABLES:
%  Problem           = structure containing optimization problem data
%    .File           =   structure of file names
%      .F            =     name of functions file
%      .O            =     name of Omega file
%      .N            =     name of discrete neighbor file
%      .C            =     name of pre-existing Cache file
%      .P            =     name of parameter file
%    .nameCache      =   name of the base workspace Cache variable
%    .fType          =   type of functions file (C=C, F=FORTRAN, M=MATLAB)
%    .Omega          =   structure defining linear constraints l <= Ax <= u
%    .iterate0       =   initial iterates
%    .isMVP          =   flag indicating an MVP problem
%    .adaptOmega     =   flag for updating an MVP's Omega parameters
%    .Param          =   structure of user-defined parameters
%      .maxNp        =   maximum number of categorical variables
%      .maxNx        =   maximum number of continuous variables
%      .maxNc        =   maximum number of nonlinear constraints
%      .pollBasis    =   basis used to construct custom Poll directions
%  Options           = structure of MADS parameters
%    .tolBind        =   tolerance for flagging active linear constraints
%    .Poll.type      =   string code indicating choice of Poll strategy
%    .Poll.order     =   string code indicating choice of Polling order
%    .Poll.center    =   number of the filter point that is Poll center
%    .delta0         =   initial poll size
%    .deltaMax       =   maximum poll size
%    .meshRefine     =   mesh refinement factor
%    .meshCoarsen    =   mesh coarsening factor
%    .hmin           =   minimum constraint violation of an infeasible point
%    .hmax           =   maximum allowable constraint violation
%    .Term           =   criteria for terminating the MADS function
%      .delta        =     lowest allowed poll size
%      .nIter        =     maximum number of iterations
%      .nFunc        =     maximum number of function calls
%      .time         =     maximum allowable CPU time
%      .nFails       =     maximum number of consecutive Poll failures
%    .EPoll.fTrigger =   objective function extended poll trigger
%    .EPoll.hTrigger =   constraint violation function extended poll trigger
%    .nSearches      =   number of different SEARCH types to be used
%    .Search(n)      =   structure containing parameters for each SEARCH
%      .type         =     string code indicating choice of SEARCH strategy
%      .label        =     long text label of SEARCH strategy
%      .nIter        =     number of iterations to perform the SEARCH
%      .nPoints      =     maximum number of points evaluated in SEARCH step
%      .complete     =     flag indicating all Search points to be evaluated
%      .sfile        =     names of optional user surrogate files
%      .file         =     name of optional file containing the user SEARCH
%    .computeGrad    =   flag for computing any available gradient
%    .accelerate     =   flag for accelerating mesh refinement
%    .useFilter      =   turns on/off the filter for nonlinear constraints
%    .runStochastic  =   flag to run as a stochastic optimization problem
%    .fixCategorical =   flag to hold categorical variable values as constant
%    .plotFilter     =   turns on/off real-time filter plot
%    .fplothandle    =   handle of the filter plot
%    .loadCache      =   flag for loading Cache of iterates
%  RunData           = structure of MADS run statistics and parameters
%    .porder         =   fixed Poll order
%    .scale          =   scale factors for scaling mesh directions
%    .stopRun        =   flag for stopping MADS run immediately
%  iterate0          = initial iterates
%  Cache             = structure of all previously evaluated iterates
%    .iterate        =   vector of all previously computed iterates
%    .Filter         =   structure of filter data
%      .plot         =     turns on/off filter plot
%      .plothandle   =     handle for the axes of the plot
%    .size           =   number of points in the Cache
%  Defaults          = defaults for all MADS parameters
%    .Options        =   Default values for Options
%    .Types          =   list of strings corresponding to possible choices
%      .poll         =     possible poll strategies
%      .pollOrder    =   possible poll order strategies
%  n0                = number of initial iterates
%  ind               = index that identifies which error was flagged
%  TempFilter        = used to reconstruct a filter from a previous run
%    .plothandle     = temporary storage for Cache.Filter.plothandle
%  Param             = user-defined parameters
%    .maxNp          =   maximum number of categorical variables
%    .maxNx          =   maximum number of continuous variables
%    .pollBasis      =   basis used to construct custom Poll directions
%    .pollOrder      =   fixed Poll order
%  iterate           = initial iterate after being processed
%*******************************************************************************
function [Problem,Options,RunData] = processInput(Problem,iterate0,varargin)

% Set Options, depending on the number of input arguments
narginchk(2,4);
Defaults = mads_defaults(Problem.typeProblem);
switch nargin
   case {2}
      Options = Defaults.Options;
   case {3,4}
      Options = deal(varargin{1});
      Options.tolBind = Defaults.Options.tolBind;
end

% Delete any abort file that exists
if exist(Options.abortFile,'file') && strcmpi(Problem.typeProblem,'TRUTH')
   delete(Options.abortFile);
end

% Set up diary file if a debugging option is selected
if Options.debug && strcmpi(Problem.typeProblem,'TRUTH')
   if exist(Problem.File.D,'file')
      diary off
      delete(Problem.File.D);
   end
   diary(Problem.File.D);
end
if Options.debug >= 1
   disp([Problem.typeProblem,' (  0):  Begin Initial Points']);
end

% Process initial points into the proper format
n0 = length(iterate0);
if isstruct(iterate0(1))
   [Problem.iterate0(1:n0).x] = deal(iterate0.x);
   [Problem.iterate0(1:n0).p] = deal(iterate0.p);
else
   [Problem.iterate0(1:n0).x] = deal(iterate0);
   [Problem.iterate0(1:n0).p] = deal({});
end

for k = 1:n0
   Problem.iterate0(k).n = length(iterate0(k).x);
end

% Option for fixing categorical variables
if Options.fixCategorical
   Problem.Param.p = Problem.iterate0(1).p;
   [Problem.iterate0(:).p] = deal({});
end

% Check input data for errors
Options.computeGrad = strncmpi(Options.Poll.type,'Gradient',8);
Problem.isMVP       = ~isempty(Problem.iterate0(1).p);
if Problem.isMVP && ~exist(Problem.File.N,'file')
   Problem.File.N = 'default_N';
   warning('mads:input','Discrete neighbors file not found.  Default is used');
end
if ~exist(Problem.File.X,'file'), Problem.File.X = []; end
checkErrors(Problem,Options,Defaults.Types);

% Adjust DACE/NW initial guess, as appropriate
for k = 1:Options.nSearches
   switch Options.Search(k).type
   case {'NW','SPS-NW'}
      if Options.nw(k).sigma < Options.nw(k).lower || ...
         Options.nw(k).sigma > Options.nw(k).upper
         Options.nw(k).sigma = mean([Options.nw(k).lower, ...
                                     Options.nw(k).upper]);
      end
   end
end

% Add a bogus empty Search to the end of the list, for convenience
if Options.nSearches > 0
   Options.Search(end+1) = struct('type','None','label','None', ...
                                  'nIter',Inf,'nPoints',0, ...
                                  'complete',1,'sfile','','file','',...
                                  'local',0,'merit',0,'param',0, ...
                                  'cblgs',Options.Search(end).cblgs);
   Options.nSearches = Options.nSearches + 1;
end

% Adjust choices for derivative options when they are unavailable
if Options.computeGrad && (Problem.fType == 'M') ...
                       && (abs(nargout(Problem.File.F)) <= 2)
   Options.computeGrad  = 0;
   Options.Poll.type    = 'Standard_2n';
   warning('mads:input', ...
           'Derivatives unavailable. Using standard 2n directions');
end

% Fix certain choices for MADS to enforce convergence theory
if ~isempty(strfind(Options.Poll.type,'MADS'))
   Options.meshCoarsen = 2;
   Options.meshRefine  = 1/2;
   Options.accelerate  = 0;
end
Options.delta0 = min(Options.delta0,Options.deltaMax);

% Change options if a filter is not used
if ~Options.useFilter
   Options.Poll.center = 0;
   Options.plotFilter = 0;
   for k = 1:Options.nSearches
      if Options.Search(k).type(2) == 'P'
         Options.Search(k).type = 'None';
      end
   end
end

% Convert Problem.File.F from cell array to string if MultiFi not selected
if ~Options.useMultiFi && iscell(Problem.File.F)
   Problem.File.F = Problem.File.F{end};
end

% Clear any previous appdata
if isappdata(0,Problem.nameCache)
   rmappdata(0,Problem.nameCache);
end
if isappdata(0,'SUR') && ~strcmp(Problem.nameCache,'sCACHE')
   rmappdata(0,'SUR');
end
if isappdata(0,'ALTCONE')
   rmappdata(0,'ALTCONE');
end

% Delete the existing history file
if exist(Problem.File.H,'file')
   delete(Problem.File.H);
end

% Initialize Cache or load a pre-existing one
Cache = createCache(Options,Problem.File.C);
setappdata(0,Problem.nameCache,Cache);

RunData.tic = cputime;

% Set memory allocation variables and user-defined parameters
Problem.maxNp = 2*length(Problem.iterate0(1).p);
Problem.maxNx = (Problem.iterate0(1).n)*(2*Problem.isMVP + ~Problem.isMVP);
if isfield(Problem,'Param')
   ParamField = {'maxNp','maxNx','PollBasis','PollOrder'};
   for k = 1:length(ParamField)
      if isfield(Problem.Param,ParamField{k})
         Problem.(ParamField{k}) = Problem.Param.(ParamField{k});
      end
   end
end

% Flag Custom Poll Order or Directions error
if strncmp(Options.Poll.type,'Custom',6) && ~isfield(Problem,'PollBasis')
   error('mads:user','User-provided Poll Basis not found (processInput).');
end
if strcmp(Options.Poll.order,'Custom') && ~isfield(RunData,'porder')
   error('mads:user','User-provided Poll Order not found (processInput).');
end

% Initialize Omega and process initial iterates
[Problem.Omega,RunData.scale] = updateOmega(Problem,Options, ...
                                            Problem.iterate0(1));
if exist(Problem.File.O,'file')
   Problem.adaptOmega = nargin(Problem.File.O) > 1;
else
   Problem.adaptOmega = 0;
end

% Process initial iterates
rng(9);
[iterate,success,TempFilter] = ...
         evalPointSet('0',Problem,Problem.iterate0,1,Options,Cache.Filter,2);
Cache = getappdata(0,Problem.nameCache);
Cache.Filter = TempFilter;
setappdata(0,Problem.nameCache,Cache)

% Create a valid Cache point if none exists
if success < 0 || all(~isfinite([Cache.iterate(1:Cache.size).f]))
   warning('mads:x0','Given initial point(s) not feasible')
   n0 = length(iterate0);
   n  = length(iterate0(1).x);
   iterate1 = dealIterate(iterate0(1),n0,zeros(n,n0));
   [iterate1.n] = deal(n);
   err = zeros(1,n0);
   for i = 1:n0
      [iterate1(i),err(i)] = makeFeasible(Problem.iterate0(i),Problem.Omega);
      if err(i) == -9
         error('mads:missingToolbox', ...
               'Optimization Toolbox required (makeFeasible)');
      end
   end
   if all(err <= 0)
      error('mads:x0','MADS cannot find feasible initial point (makeFeasible)');
   end
   [iterate,~,TempFilter] = ...
        evalPointSet('0',Problem,iterate1,1,Options,Cache.Filter,2); 
   Cache = getappdata(0,Problem.nameCache);
   Cache.Filter = TempFilter;
   setappdata(0,Problem.nameCache,Cache);
end
Problem.maxNc = length(iterate(1).c)*(2*Problem.isMVP + ~Problem.isMVP);

% Initialize RunData parameters
if nargin >= 4                      % For resuming previous run
   RunData = varargin{2};
   RunData.stopRun = 0;
else
   [Problem,RunData] = update(Problem,Options,RunData,-9);
end
if Options.debug >= 1
   disp([Problem.typeProblem,' (  0):  End Initial Points']);
end
return

%*******************************************************************************
% checkErrors:  Check for errors in the input data.
% ------------------------------------------------------------------------------
% Called by: processInput
% VARIABLES:
%  Problem           = structure containing optimization problem data
%    .File           =   structure of file names
%      .F            =     name of functions file
%      .O            =     name of Omega file
%      .N            =     name of discrete neighbor file
%      .C            =     name of pre-existing Cache file
%      .P            =     name of parameter file
%    .nameCache      =   name of the base workspace Cache variable
%    .fType          =   type of functions file (C=C, F=FORTRAN, M=MATLAB)
%    .Omega          =   structure defining linear constraints l <=Ax <= u
%    .iterate0       =   initial iterates
%    .isMVP          =   flag indicating an MVP problem
%    .adaptOmega     =   flag for updating an MVP's Omega parameters
%    .maxNp          =   maximum number of categorical variables
%    .maxNx          =   maximum number of continuous variables
%    .maxNc          =   maximum number of nonlinear constraints
%    .pollBasis      =   basis used to construct custom Poll directions
%  Options           = structure of MADS parameters
%    .tolBind        =   tolerance for flagging active linear constraints
%    .pollStrategy   =   string code indicating choice of Poll strategy
%    .pollOrder      =   string code indicating choice of Polling order
%    .pollCenter     =   number of the filter point that is Poll center
%    .delta0         =   initial poll size
%    .deltaMax       =   maximum poll size
%    .meshRefine     =   mesh refinement factor
%    .meshCoarsen    =   mesh coarsening factor
%    .hmin           =   minimum constraint violation of an infeasible point
%    .hmax           =   maximum allowable constraint violation
%    .Term           =   criteria for terminating the MADS function
%      .delta        =     lowest allowed poll size
%      .nIter        =     maximum number of iterations
%      .nFunc        =     maximum number of function calls
%      .time         =     maximum allowable CPU time
%      .nFails       =     maximum number of consecutive Poll failures
%    .EPoll.fTrigger =   objective function extended poll trigger
%    .EPoll.hTrigger =   constraint violation function extended poll trigger
%    .nSearches      =   number of different SEARCH types to be used
%    .Search(n)      =   structure containing parameters for each SEARCH
%      .type         =     string code indicating choice of SEARCH strategy
%      .label        =     long text label of SEARCH strategy
%      .nIter        =     number of iterations to perform the SEARCH
%      .nPoints      =     maximum number of points evaluated in SEARCH step
%      .sfile        =     name of optional user surrogate files
%      .file         =     name of optional file containing the user SEARCH
%    .computeGrad    =   flag for computing any available gradient
%    .accelerate     =   flag for accelerating mesh refinement
%    .useFilter      =   turns on/off the filter for nonlinear constraints
%    .plotFilter     =   turns on/off real-time filter plot
%    .fplothandle    =   handle of the filter plot
%    .loadCache      =   flag for loading Cache of iterates 
%  indNone           = indices of Searches that are of type "None" 
%  errorMsg          = cell array of error messages
%  errorCheck        = vector of flags that track input data errors
%*******************************************************************************
function checkErrors(Problem,Options,Types)

% Check input parameters for errors and display, if any
errorMsg{1}    = 'Required Functions file not found';
errorMsg{2}    = 'Required Omega file not found';
errorMsg{3}    = 'Invalid choice of Search Method';
errorMsg{4}    = 'Invalid choice of Poll Directions';
errorMsg{5}    = 'Invalid choice of Polling Order';
errorMsg{6}    = 'Invalid choice of Poll Center';
errorMsg{7}    = 'Termination criteria must be positive numbers or infinity';
errorMsg{8}    = 'Termination criteria cannot all be infinite';
errorMsg{9}    = 'Initial Mesh Size must be a strictly positive number';
errorMsg{10}   = 'Maximum Mesh Size must be greater than Initial Mesh Size';
errorMsg{11}   = 'Mesh Refinement Factor must be a rational number in (0,1)';
errorMsg{12}   = 'Mesh Coarsening Factor must be a rational number > 1';
errorMsg{13}   = 'Min Filter Constraint Violation must be nonnegative';
errorMsg{14}   = 'Max Filter Constraint Violation must be positive';
errorMsg{15}   = 'Filter Constraint Violation Values must have Min < Max';
errorMsg{16}   = 'Extended Poll Triggers must be positive numbers';
errorMsg{17}   = 'Number of runs must be positive';
errorMsg{18}   = 'Lower bound on surrogate parameter exceeds its upper bound';
errorMsg{19}   = 'A Nadaraya-Watson surrogate requires a stochastic problem.';
errorCheck(1)  = ~exist([Problem.File.F, '.', lower(Problem.fType)],'file');
errorCheck(2)  =  exist(Problem.File.N,'file') && ~exist(Problem.File.O,'file');
errorCheck(3)  = Options.nSearches > 0 && ...
                 any(~ismember({Options.Search(:).type}, Types.search));
errorCheck(4)  = ~ismember(Options.Poll.type,  Types.poll);
errorCheck(5)  = ~ismember(Options.Poll.order, Types.pollOrder);
errorCheck(6)  = ~isfinite(Options.Poll.center) || Options.Poll.center < 0;
Term           = [Options.Term.delta,Options.Term.nIter, ...
                  Options.Term.nFunc,Options.Term.time];
errorCheck(7)  = any(isnan(Term)) || any(Term < 0);
errorCheck(8)  = all(isinf(Term));
errorCheck(9)  = ~isfinite(Options.delta0) || Options.delta0 <= 0;
errorCheck(10) = Options.delta0 > Options.deltaMax;
errorCheck(11) = ~isfinite(Options.meshRefine) ...
              || Options.meshRefine <= 0 || Options.meshRefine >= 1;
errorCheck(12) = ~isfinite(Options.meshCoarsen) || Options.meshCoarsen < 1;
errorCheck(13) = ~isfinite(Options.hmin)        || Options.hmin <  0;
errorCheck(14) = ~isfinite(Options.hmax)        || Options.hmax <= 0;
errorCheck(15) = Options.hmax <= Options.hmin;
errorCheck(16) = isnan(Options.EPoll.fTrigger) || Options.EPoll.fTrigger <=0 ...
              || isnan(Options.EPoll.hTrigger) || Options.EPoll.hTrigger <=0;
errorCheck(17) = Options.runStochastic && ~isfinite(Options.Stoch.nRuns) <= 0;
for k = 1:Options.nSearches
   switch Options.Search(k).type
      case {'NW','SPS-NW'}
         errorCheck(18) = Options.nw(k).lower > Options.nw(k).upper;
         errorCheck(19) = ~Options.runStochastic;
      otherwise
         errorCheck(18) = 0;
   end
   if errorCheck(18), break, end
end

ind = find(errorCheck);
if ~isempty(ind)
   error('mads:input',[errorMsg{ind(1)},' (checkErrors).']);
end
return

%*******************************************************************************
% createCache:  Create and initialize the Cache.
% ------------------------------------------------------------------------------
% Called by: processInput
% Calls:     createFilter
% VARIABLES:
%  Cache            = structure containing a Cache of iterates
%    .tol           =   tolerance for an iterate being in the cache
%    .iterate       =   vector of iterates
%      .n           =     size of the continuous variable space
%      .x(nc)       =     values of continuous variables 
%      .p(nd)       =     values of categorical variables
%      .f           =     objective function value (f-value)
%      .c(m)        =     constraint function values (c(x) <= 0)
%      .h           =     constraint violation function value (h-value)
%      .gradf(nc)   =     partial derivatives of f (if available)
%      .gradc(m,nc) =     partial derivatives of c (if available)
%      .gradh       =     partial derivatives of h (if available)
%      .param       =     unspecified user-determined values
%    .size          =   number of iterates
%    .xnorm         =   vector of inf-norms of the iterates
%    .pID           =   vector of unique IDs for categorical variables
%    .isSurrPoint   =   vector of flags of points used in surrogate recal
%    .bfp           =   vector of best cumulative feasible point f-values
%    .lip           =   matrix of least infeasible point f- and h-values
%    .nHits         =   number of Cache hits
%    .Filter        =   structure containing filter data
%  Options          = structure containing SEARCH and POLL parameters
%    .computeGrad   =   logical indicating gradients are to be computed
%    .hmin          =   minimum constraint violation of an infeasible point
%    .hmax          =   maximum allowable constraint violation
%    .plotFilter    =   logical for plotting real-time filter
%    .fplothandle   =   handle for the real-time filter plot
%    .tolCache      =   tolerance for a point to be a Cache point
%    .CacheFile     =   name of file containing a pre-existing Cache
%*******************************************************************************
function Cache = createCache(Options,CacheFile)

% Create a new temporary, unpopulated Cache -- it becomes Cache if none exists
TempCache     = [];
TempCache.tol = Options.tolCache;
TempCache.iterate(1).x     = [];
TempCache.iterate(1).p     = {};
TempCache.iterate(1).n     = 0;
TempCache.iterate(1).type  = ' ';
TempCache.iterate(1).sf    = [];
TempCache.iterate(1).sc    = [];
TempCache.iterate(1).sh    = [];
TempCache.iterate(1).f     = [];
TempCache.iterate(1).c     = [];
TempCache.iterate(1).h     = [];
TempCache.iterate(1).param = [];
if Options.computeGrad
   TempCache.iterate(1).sgradf = [];
   TempCache.iterate(1).sgradc = [];
   TempCache.iterate(1).sgradh = [];
   TempCache.iterate(1).gradf  = [];
   TempCache.iterate(1).gradc  = [];
   TempCache.iterate(1).gradh  = [];
end
TempCache.size                = 0;
TempCache.xnorm               = [];
TempCache.pID                 = [];
TempCache.isSurrPoint         = [];
TempCache.nFunc               = [];
TempCache.sdev                = [];
TempCache.bfp                 = [];
TempCache.lip                 = [];
TempCache.nHits               = 0;
TempCache.iterate(1)          = [];            % This is a key statement!
TempCache.Filter = createFilter(Options.hmin,Options.hmax,0,...
                                Options.plotFilter,Options.fplothandle);

% A Cache file exists and will be used
if Options.loadCache && exist(CacheFile,'file')
   load(CacheFile,'Cache');

   % Ensure the existing Cache has the proper fields
   try
      for fname = fieldnames(TempCache)
         if ~isfield(Cache,fname)
            Cache.(fname) = TempCache.(fname);
         end
      end
      Cache = orderfields(Cache,TempCache);
   catch
   end

   % Setup existing filter
   TempFilter = createFilter(0.1,1,0,1,Options.fplothandle);
   Cache.Filter.plothandle = TempFilter.plothandle;
   if Cache.Filter.plot
      plotFilter(Cache.Filter,Cache.Filter.hmax,Cache.iterate,10);
   end
else
   Cache = TempCache;
end
return

%*******************************************************************************
% createFilter:  Create and initialize a filter.
% ------------------------------------------------------------------------------
% Called by: processInput, mvpPoll
% VARIABLES:
%  Filter        = structure containing the filter
%    .feasible   =   indices of feasible iterates, sorted by f-value
%    .F          =   indices of iterates in the filter, sorted by h-value
%    .hmin       =   minimum allowable h-value to be in Filter.F
%    .hmax       =   maximum allowable h-value to be in Filter.F
%    .strict     =   flag indicating that hmax is strictly enforced
%    .plot       =   flag for plotting real-time filter
%    .plothandle =   handle of the real-time filter plot
%  hmin          =   input value for Filter.hmin
%  hmax          =   input value for Filter.hmax
%  strict        =   input value for Filter.strict
%  plot          =   input value for Filter.plot
%  plothandle    =   input value for Filter.plothandle
%*******************************************************************************
function Filter   = createFilter(hmin,hmax,strict,fplot,plothandle)
if hmin >= hmax
   error('mads:filter',['Attempt to create filter with ',...
                        'incompatible hmin/hmax (createFilter).']);
end
Filter.hmin       = hmin;
Filter.hmax       = hmax;
Filter.strict     = strict;
Filter.plot       = fplot;
Filter.F          = [];
Filter.feasible   = [];
Filter.BFP        = [];
Filter.LIP        = [];
Filter.IIP        = [];
if isempty(plothandle) && fplot
   figure; plothandle = gca;
   set([plothandle; get(plothandle,'Children')],'Visible','off');
end
Filter.plothandle = plothandle;
return

%*******************************************************************************
% makeFeasible:  Projects an infeasible iterate (w.r.t. Omega) into Omega.
% ------------------------------------------------------------------------------
% Called by: processInput
% Calls:     inOmega, runLP
% VARIABLES:
%  iterate    = current iterate that will be made feasible
%    .x       =   continuous variables
%    .p       =   categorical variables
%  err        = error code indicating success/failure in finding feasible point
%  Omega      = structure defining linear constraints, l <= A*x <= u
%    .A       =   matrix of linear coefficients
%    .l       =   vector of lower bounds
%    .u       =   vector of upper bounds
%    .plist   =   list of allowable discrete variable values
%  pass       = flag indicating if categorical variables are feasible
%  [m,n]      = dimensions of linear constraint matrix Omega.A
%  ind        = indices of infeasible lower and upper bounds
%  c,A,b      = matrices that define the LP, min {c'x: Ax <= b}
%  OptOptions = structure used in LP solver to set optimization options
%  dist       = closest distance of iterate to Omega
%  err        = error flag for LP solver
%*******************************************************************************
function [iterate,err] = makeFeasible(iterate,Omega)

% Adjust categorical variables
for k = 1:length(iterate.p)
   if ischar(iterate.p{k})
      pass = any(strcmp(iterate.p{k},Omega.plist{k}));
   else
      pass = any([Omega.plist{k}{:}] == iterate.p{k});
   end
   if ~pass
      iterate.p{k} = Omega.plist{k}{1};
   end
end

% Adjust bounds if necessary
[m,n] = size(Omega.A);
if n <= m && isequal(Omega.A(1:n,:),eye(n))
   ind = find(iterate.x - Omega.l(1:n) < 0.0);
   iterate.x(ind) = deal(Omega.l(ind));
   ind = find(iterate.x - Omega.u(1:n) > 0.0);
   iterate.x(ind) = deal(Omega.u(ind));
   if inOmega(iterate,Omega)
      err = 1;
      return
   end
end

% Set up linear program and remove non-finite rows
c = [zeros(n,1); 1];
b = [iterate.x; -iterate.x; Omega.u; -Omega.l];
A = [ eye(n) ,  -ones(n,1); ...
     -eye(n) ,  -ones(n,1); ...
      Omega.A,  zeros(m,1); ...
     -Omega.A,  zeros(m,1)];
ind      = find(~isfinite(b));
b(ind)   = [];
A(ind,:) = [];

% Solve LP and assess if new iterate is feasible
[iterate.x,~,err] = runLP(c,A,b);
iterate.x(end) = [];
return

%*******************************************************************************
% FUNCTIONS CALLED BY MADS MAIN LOOP
%*******************************************************************************
%*******************************************************************************
% search:  Performs the SEARCH step of the MADS/GPS algorithm.
%          A finite number of points are evaluated according to the k-th Search
%          strategy set by the variable Search(k).Type. (Default: no search)
% ------------------------------------------------------------------------------
% Called by: mads
% Calls:     lhs, ga, recalSurrogate, optimizeSurrogate, evalPointSet, poll, 
%            updateOmega, gridsamp (DACE)
% VARIABLES:
%  success         = flag indicating success or failure of Poll step
%  Filter          = structure of filter data
%    .F            =   indices of Cache points in the filter
%  Problem         = structure of optimization problem data
%    .nameCache    =   name of the base workspace Cache variable
%    .Omega        =   structure defining linear constraints
%      .plist      =     list of allowable discrete variable values
%    .isMVP        =   flag indicating an MVP problem
%    .adaptOmega   =   flag for updating an MVP's Omega parameters
%  Options         = structure of MADS parameters
%    .pollComplete =   turns on/off evaluation of ALL Poll points
%    .computeGrad  =   flag for computing any available gradients
%    .SurOptimizer =   string containing name of surrogate optimizer
%    .Search       =   vector of structures of Search data
%    .nSearches    =   number fo different Searches used during run
%    .pollStrategy =   string identifying selected Poll strategy
%    .pollComplete =   flag for performing a complete Poll step
%  RunData         = structure of run statistics and parameters
%    .scale        =   scale factors for scaling of mesh directions
%    .pcenter      =   current Poll center
%      .n          =     number of continuous variables
%    .delta        =   poll size parameter
%    .nIter        =   iteration count
%  Cache           = structure containing of data on all past iterates
%    .iterate      =   vector of previously processed iterates
%    .size         =   number of iterates in Cache
%    .tol          =   tolerance for determining if a point is in Cache
%  full            = temporary storage of Options.Poll.complete
%  grad            = temporary storage of Options.computeGrad
%  optimizer       = temporary storage of Options.SurOptimizer
%  Cache           = temporary storage of Cache.iterate
%  tol             = temporary storage of Cache.tol
%  Search          = structure of Search data 
%    .type         =   string indicating which search method to use
%    .nIter        =   number of iterations to perform Search
%    .nPoints      =   number of search points to be evaluated
%    .sfile        =   names of optional user surrogate files
%    .file         =   name  of optional user-defined Search file
%    .dace         =   structure of DACE surrogate data
%      .reg        =     regression function handle
%      .corr       =     correlation function handle
%      .theta      =     correlation function parameter
%      .isotropic  =     flag indicating isotropic (all equal) theta values
%    .nw           =   structure of NW surrogate data
%      .kernel     =     name of NW kernel function
%      .local      =     flag indicating if surrogate is restricted to a region
%  surrogate       =   structure of surrogate data
%    .recalibrator =     surrogate recalibration function filename
%    .evaluatror   =     surrogate evaluation function filename
%  pcenter         = temporary poll center
%  D               = Mesh directions used in LHS
%  S               = Search points to be evaluated
%  N               = discrete neighbors of LHS points
%  q               = number of mesh points in each direction
%  LB,UB           = lower/upper bounds on continuous variables
%  nS              = number of search points to be evaluated
%  meshPoints      = mesh points for coarse mesh point evaluations
%  param           = cell array of parameters output from user surrogate file
%*******************************************************************************
function [success,Filter,term] = search(Problem,Options,RunData,Filter)

isSurrogate = @(x) any(endsWith(upper(x),["DACE","RBF","NW","CUSTOMS"]));
Cache = getappdata(0,Problem.nameCache);
success = -1;
term = 0;

% Shortcuts
grad    = Options.computeGrad;
tol     = Cache.tol;
ind     = 1:Cache.size;
pcenter = RunData.pcenter;
n       = pcenter.n;
LB      = Problem.Omega.l(1:n);
UB      = Problem.Omega.u(1:n);
nIter   = RunData.nIter;
debug   = Options.debug;
msg     = ' requires finite lower and upper bounds on every variable';

% Identify which Search will be performed
if isempty(Options.Search), return, end
Search           = Options.Search(1);
Search.dace      = Options.dace(1);
Search.nw        = Options.nw(1);
Search.rbf       = Options.rbf(1);
Search.cmaes     = Options.cmaes(1);
Search.optimizer = Options.SurOptimizer;
for k = 2:Options.nSearches
   if RunData.nIter >= sum([Options.Search(1:k-1).nIter])
      Search = Options.Search(k);
      if k < Options.nSearches
         Search.dace  = Options.dace(k);
         Search.nw    = Options.nw(k);
         Search.rbf   = Options.rbf(k);
         Search.cmaes = Options.cmaes(k);
         Search.optimizer = Options.SurOptimizer;
      end
   end
end
if strcmpi(Search.type,'NONE') || ~Search.nPoints, return, end

% Make sure any user search files exist
if any(startsWith(Search.type,{'SPS','CUS'})) && ...
      (isempty(Search.file) || ~exist(func2str(Search.file),'file'))
   error('mads:search',['User file (',func2str(Search.file), '.m) not found.']);
end

% Option to include only iterates with same categorical variable values (Dunlap)
plist = {};
if Problem.isMVP
   plist = Problem.Omega.plist;
   pUB   = cellfun(@length,plist);
   if Options.mvp1Surrogate
      ind = find([Cache.iterate.n] == pcenter.n);
   else
      if strcmp(Search.optimizer,'cmaes')
         error('mads:search','CMA-ES cannot be used with multiple surrogates')
      end
      pID = getpID(pcenter.p,plist);
      ind = find(Cache.pID == pID);
   end
   Cache.iterate     = Cache.iterate(ind);
   Cache.isSurrPoint = Cache.isSurrPoint(ind);
   ind = 1:length(Cache.iterate);
end

% Change Surrogate to LHS if Cache is empty and filter out points for surrogates
switch upper(Search.type)
   case {'DACE','NW','RBF','SPS-DACE','SPS-NW','SPS-RBF','CUSTOMS'}
      if strcmp(Search.optimizer,'cmaes') && length(Problem.Omega.u) < n
         error('mads:search',['CMA-ES software',msg])
      end
      ind = find(Cache.isSurrPoint);
      if length(Cache.iterate) <= 1
         if isSurrogate(Options.Search(1).type)
            Search = struct('type','LHS','label','Latin Hypercube Search', ...
                            'nIter',1,'nPoints',Options.lhsStrength*pcenter.n, ...
                            'complete',1,'sfile',{'',''},'file','');
         else
            Search = Options.Search(1);
         end
      end

   % Retrieve linear constraints including variable bounds
   case {'CCD','LHS','MESH','CMAES','PS'}
      Options.removeRedundancy = 0;
      [Problem.Omega,RunData.scale] = updateOmega(Problem,Options,pcenter);
      LB = Problem.Omega.l(1:n);
      UB = Problem.Omega.u(1:n);
      I  = eye(n);
end
Cache = Cache.iterate(ind);

% Set sampling bounds (same as problem bounds unless provided by user)
if isfield(Problem,'Param') && isfield(Problem.Param,'dx0')
   sLB = Problem.Param.iterate0.x - Problem.Param.dx0;
   sUB = Problem.Param.iterate0.x + Problem.Param.dx0;
else
   sLB = LB;
   sUB = UB;
end

% Call appropriate Search functions
if Options.debug >= 2
   disp(['   Search type:  ', Search.type]);
   disp(['   # Points:     ', int2str(Search.nPoints)]);
end
switch upper(Search.type)

% Poll around best n filter points
case {'SPOLLI','GPOLLI'}
   if startsWith(Search.type,'G') && Options.computeGrad
      Options.Poll.type = 'Gradient_3n_L2';
   else
      Options.Poll.type = 'Standard_2n';
   end
   Options.Poll.complete = Search.complete;
   for k = 1:min(length(Filter.F),Search.nPoints)
      pcenter = Cache(Filter.F(k));
      if Problem.adaptOmega
         [Problem.Omega,RunData.scale] = updateOmega(Problem,Options,pcenter);
      end
      [success,Filter] = poll('S',Problem,Options,RunData,pcenter,Filter);
      if success > -1 && ~Options.Poll.complete, break, end
   end
   return

% Central Composite Design
case 'CCD'

   % Error if any variable has non-finite bounds
   if ~isequal(Problem.Omega.A(1:n,1:n),I) || ~all(isfinite(sLB.*sUB))
      error('mads:search',['CCD',msg]);
   end

   % Form generic CCD
   Z = [-1; 1];
   for k = 1:n-1
      e = ones(size(Z,1),1);
      Z = [Z, -e; Z, e];  %#ok
   end
   Z = [zeros(1,n); sqrt(n)*I; -sqrt(n)*I; Z];

   % Translate CCD to the sampling bound constrained region
   nS = size(Z,1);
   [S(1:nS).x] = deal(pcenter.x);
   [S(1:nS).p] = deal(pcenter.p);
   for j = 1:nS
      S(j).x = sLB + ((sUB - sLB)/(2*sqrt(n))).*(Z(j,:)' + sqrt(n));
   end

% Latin Hypercube Search on the mesh
case 'LHS'

   % Error if any variable has non-finite bounds
   if ~isequal(Problem.Omega.A(1:n,1:n),I) || ~all(isfinite(sLB.*sUB))
      error('mads:search',['LHS',msg]);
   end

   % Augment with categorical variables
   if Problem.isMVP
      sLB  = [sLB; ones(length(Problem.Omega.plist),1)];
      sUB  = [sUB; pUB(:)];
   end

   % Loop for generating enough valid LHS points
   nPoints  = Search.nPoints;
   S        = [];
   nS       = 0;
   addNoise = 0;
   nFail    = 10;
   iFail    = 0;
   if nPoints == 1
      nPoints = Options.lhsStrength*(n+length(pcenter.p));
   end
   while nS < nPoints

      % Generate LHS points in the feasible region
      nGrid = nPoints*Options.lhsBinFactor;
      z     = lhs(sLB,sUB,nPoints,nGrid,addNoise,plist);

      % Screen LHS points for closed constraint feasibility
      S1 = dealIterate(pcenter,nPoints,z,plist);
      S1 = screenClosedConstraints(S1,Problem.File.X,nPoints);

      % Add LHS points to search points or throw error if no valid points found
      if isempty(S1)
         iFail = iFail + 1;
         if iFail > nFail
            error('mads:search:LHS','No X-feasible points found by LHS');
         end
      else
         S  = [S, S1];  %#ok
         nS = length(S);
         if Options.debug >= 2
            disp(['   Points/Valid: ',int2str(nPoints),'/',int2str(nS)]);
         end
         iFail = 0;
      end
   end

% Deterministic Sampling on a Coarse Mesh (requires DACE package)
case 'MESH'

   % Error if any variable has non-finite bounds
   if ~isequal(Problem.Omega.A(1:n,1:n),I) || ~all(isfinite(sLB.*sUB))
      error('mads:search',['MESH',msg]);
   end

   % Call the DACE mesh function to construct mesh points
   nGrid = ceil(Search.nPoints^(1/n));
   if Problem.isMVP
      nGrid = [nGrid*ones(1,pcenter.n), pUB(:)'];
   end
   meshPoints = gridsamp([sLB,sUB]',nGrid);
   S = dealIterate(pcenter,size(meshPoints,1),meshPoints');

case 'PS'
   Options.vTol = Options.Term.delta;
   Options.maxIterations = 1;
   Options.nPoints = Search.nPoints;
   [x,fx,term] = ps(Problem.File.F,[LB,UB],pcenter.p,Options);
   nS = length(fx);
   S  = dealIterate(pcenter,nS,x);
   [S(1:nS).f] = deal(fx);

% Covariance Matrix Adaptive Evolutionary Strategy
case 'CMAES'
   penalty = setPenaltyTerm(Cache,Filter);
   if exist(Problem.File.O,'file')
      [A,lb,ub] = feval(Problem.File.O,n);
   end
   f      = @(x) eaPenalty(x,pcenter.p,Problem.File.F,A,lb,ub,penalty);
   [X,fX] = ea(f,LB,UB,Options,Search);
   nS     = length(fX);
   S      = dealIterate(pcenter,nS,X);
   [S(1:nS).f] = deal(fX);

% DACE Surrogate (used for both traditional and simplified physics)
case {'DACE','SPS-DACE'}
   Sites     = getDataSites(Problem,Cache,Search,Options.mvp1Surrogate, ...
                            Options.digits,debug);
   surrogate = daceSurrogate(Search,Sites,Filter,Problem.Omega,nIter,debug);
   [S,term]  = optimizeSurrogate(Problem,pcenter,surrogate,Options,Search,tol);
   
% Radial basis function or Nadaraya-Watson surrogate
case {'NW','SPS-NW','RBF','SPS-RBF','CUSTOMS'}
   Sites     = getDataSites(Problem,Cache,Search,Options.mvp1Surrogate, ...
                            Options.digits,debug);
   surrogate = constructSurrogate(Search,Problem.Omega,nIter);
   surrogate = updateSurrogateData(Sites,surrogate,Filter);
   for k = 1:length(surrogate.resp)
      surrogate.(surrogate.resp{k}) = feval(surrogate.recalibrator, ...
                                            surrogate.X,surrogate.Y(:,k), ...
                                            surrogate.arg{:});
   end
   setappdata(0,'SUR',surrogate);
   [S,term] = optimizeSurrogate(Problem,pcenter,surrogate,Options,Search,tol);

% Custom Search method
case 'CUSTOM'
   S = feval(Search.file,Problem,Options,Search,RunData.delta,pcenter);

% Invalid Search option
otherwise
   error('mads:search:choice','Invalid choice of Search Method (search).');
end

% Run CBLGS to augment with surrograte model refinement points
if isSurrogate(Search.type) && Search.cblgs.nGoal
   C = cblgs(sLB,sUB,plist,pcenter,surrogate(end),Search.cblgs,Problem.File.X);
   S = [S; C'];
end

% Evaluate the set of points generated by user-selected Search option
[~,success,Filter] = evalPointSet('S',Problem,S,Search.complete,Options,Filter,2);
return

%*******************************************************************************
% poll: Performs the POLL step of the MADS algorithm.  Given a Poll center, a
%       set of Poll direction vectors, and a mesh size parameter, the
%       neighboring mesh points are evaluated, seeking an improved mesh point.
% ------------------------------------------------------------------------------
% Called by: mads, Search, mvpPoll
% Calls:     standDirections, gradDirections, getPollOrder, evalPointSet
% VARIABLES:
%  success         = flag indicating success or failure of Poll step
%  Filter          = structure describing the filter
%  ptype           = type of poll (P=poll, E=extended poll, S=from Search step)
%  Problem         = structure describing the optimization problem
%    .Omega        =   structure defining Omega = {x: l <= A*x <= u}
%    .adaptOmega   =   flag for updating an MVP's Omega parameters
%  Options         = structure containing various MADS parameters
%    .computeGrad  =   flag for computing any available gradients
%    .pollStrategy =   string identifying selected Poll strategy
%    .pollOrder    =   string identifying selected Poll order strategy
%    .pollComplete =   flag for performing a complete Poll step
%    .tolBind      =   active constraint tolerance
%  RunData         = structure of run statistics and parameters
%    .scale        =   scale factors for scaling mesh directions
%    .delta        =   poll size parameter
%    .goodD        =   index of successful Poll direction
%    .porder       =   current poll order
%    .nFails       =   current number of consecutive Poll failures
%  pcenter         = center of Poll set, around which Poll is performed
%    .x            =   continuous variables
%    .p            =   categorical variables
%    .h            =   constraint violation function value
%    .gradf        =   gradient of f
%    .gradh        =   gradient of h
%  D               = matrix whose columns form the direction vectors
%  g               = gradient direction with which to prune
%  infGrad         = indices of non-finite partial derivatives
%  surrogate       = structure of surrogate model data
%  order           = order in which directions will be polled
%  P               = the Poll set
%  nP              = number of points in the Poll set
%*******************************************************************************
function [success,Filter,varargout] = poll(ptype,Problem,Options,RunData,...
                                           pcenter,Filter)

% A gradient poll without gradients available becomes a standard poll
if Options.computeGrad && ~isfield(pcenter,'gradf')
   Options.Poll.type = 'Standard_2n';
end
Options.tolBind = RunData.delta_m;

% Retrieve scaled Poll directions
switch Options.Poll.type
case {'Standard_2n','Standard_n+1','Custom_2n','Custom_n+1'}
   D = standDirections(pcenter.x,Problem.Omega,Options.Poll.type, ...
                       Options.tolBind,RunData.scale,...
                       Options.degeneracyScheme,~RunData.nFails);
case {'MADS_2n','MADS_n+1','OrthoMADSr_2n','OrthoMADSr_n+1','MADSr_2', ...
      'OrthoMADS_2n','OrthoMADS_n+1','MADS_2'}
   [D,RunData.bSet] = madsDirections(pcenter.x,Options.Poll.type, ...
                      Options.delta0,RunData.delta_m,RunData.scale,  ...
                      RunData.bSet);
   if Options.dimSensor
      if mod(pcenter.n,Options.dimSensor)
         error('poll:mads:badSP','Sensor dimension does not divide into x.');
      end
      D = D(:,1);
      nSensors = pcenter.n/Options.dimSensor;
      for k = 1:nSensors
         ind = 1+(k-1)*Options.dimSensor : k*Options.dimSensor;
         sensorD = madsDirections(pcenter.x(ind),Options.Poll.type, ...
                   Options.delta0,RunData.delta_m, RunData.scale(ind), ...
                   RunData.bSet);
         tempD = zeros(pcenter.n, size(sensorD,2));
         tempD(ind,:) = sensorD;
         D = [D, tempD];  %#ok
      end
   end
case {'Gradient_2n','Gradient_n+1'}
   D = standDirections(pcenter.x,Problem.Omega,Options.Poll.type, ...
                       Options.tolBind,RunData.scale, ...
                       Options.degeneracyScheme,~RunData.nFails);
   if pcenter.h <= Options.hmin
      g = pcenter.gradf;
   else
      g = pcenter.gradh;
   end
   infGrad = ~isfinite(g);
   g(infGrad) = 0;
   D(:,(infGrad'*D == 0 & g'*D > 0.0)) = [];
case {'Gradient_3n_L1','Gradient_3n_L2','Gradient_3n_LInf','Gradient_3n2n'}
   D = [];
   if any(pcenter.gradf) || (isempty(pcenter.gradh) || any(pcenter.gradh))
      D = gradDirections(pcenter,Problem.Omega,Options.Poll.type,Filter, ...
                         RunData.delta_m,Options.tolBind,RunData.scale, ...
                         Options.degeneracyScheme,RunData.success);
   end
otherwise
   error('mads:poll:choice','Invalid Poll strategy type (poll).');
end

% Construct the POLL set
nP = size(D,2);
P  = dealIterate(pcenter, nP, repmat(pcenter.x,1,nP) + RunData.delta_m*D);

% Get surrogate information for setting surrogate Poll order
surrogate = [];
if isappdata(0,'SUR')
   surrogate = getappdata(0,'SUR');
end
if strcmp(Options.Poll.order,'Surrogate') && isempty(surrogate)
   Options.Poll.order = 'DynamicRanked';
end

% Set the polling order, and evaluate the POLL points
order = getPollOrder(Options.Poll.order,D,P,RunData,surrogate);
P     = P(order);

Problem.adaptOmega = 0;
[P,success,Filter] = evalPointSet(ptype,Problem,P,Options.Poll.complete,...
                                                  Options,Filter,1);

% Set dynamic poll ordering parameters by finding which direction led to success
if success > 0
   [~,ind]  = min([P.f]);
   RunData.goodD = D(:,ind);
   RunData.porder = order;
   success = 1;
end   
varargout = {RunData};
return

%*******************************************************************************
% mvpPoll:  Performs discrete neighbor and extended polls for MVP problems.
% ------------------------------------------------------------------------------
% Called by: mads
% Calls:     evalPointSet, poll, getPollCenter, createFilter, updateOmega
% VARIABLES:
%  success          = flag indicating success or failure of MVP Poll step
%  Problem          = structure describing the optimization problem
%    .nameCache     =   name of the base workspace Cache variable
%    .File.N        =   name of discrete neighbor file
%    .Omega         =   structure defining linear constraints
%      .plist       =     list of allowable discrete variable values
%    .adaptOmega    =   flag for updating an MVP's Omega parameters
%  Options          = structure containing various MADS parameters
%    .pollComplete  =   flag for performing a complete Poll
%    .NPollComplete =   flag for performing a complete discrete neighbor Poll
%    .EPollComplete =   flag for performing a complete Extended Poll
%    .pollCenter    =   code for identifying the Poll center
%  RunData          = structure containing MADS run statistics
%    .delta         =   current poll size
%    .pcenter       =   current Poll center
%    .fxi           =   f-value threshold for triggering Extended Poll
%    .scale         =   scale factors for scaling of mesh directions
%  Cache            = database of all previously computed iterates
%    .Filter        =   structure of parameters describing the main filter
%      .hmin        =     minimum h-value of an infeasible point
%      .hmax        =     maximum allowable h-value
%  full             = temporary storage of Options.Poll.complete
%  N                = discrete neighbors of the current poll center
%    .f             =   objective function value
%    .h             =   constraint violation function value
%  BestF            = best feasible iterate
%    .f             =   objective function value
%  BestI            = least infeasible iterate
%    .h             =   constraint violation function value
%  fplusxi          = BestF.f + RunData.fxi
%  hplusxi          = BestF.h + RunData.hxi (or hmax if this is too large)
%  ePollF           = indices of N triggering Extended Poll due to f
%  ePollH           = indices of N triggering Extended Poll due to h
%  ePoll            = indices of N triggering Extended Poll due to f or h
%  new              = temporary storage
%  order            = sorted index for N by f-value or h-value
%  pcenter          = current Extended Poll center
%  nx               = number of continuous variables in pcenter
%  Filter           = structure describing the temporary local filter
%  ePsuccess        = flag indicating a successful Extended Poll step
%  unfiltered       = flag indicating an unfiltered point has been found
%*******************************************************************************
function success = mvpPoll(Problem,Options,RunData)

% Retrieve set of discrete neighbors
N = feval(Problem.File.N,Problem,RunData.pcenter,Problem.Omega.plist, ...
          RunData.delta);
if isempty(N)
   success = -1;
   return
end

% Evaluate set of discrete neighbors
Cache = getappdata(0,Problem.nameCache);
full  = Options.EPoll.completeN;
[N,success,TempFilter] = evalPointSet('N',Problem,N,full,Options,Cache.Filter,1);
Cache = getappdata(0,Problem.nameCache);
Cache.Filter = TempFilter;
setappdata(0,Problem.nameCache,Cache);
N = N(~isnan([N.f]));
%success = ~~success;

% If unsuccessful, begin EXTENDED POLL step around "good" discrete neighbors
if success < 0

   % Determine candidate poll centers for extended polling
   switch Options.useFilter
      case 0
         centerCode = 0;
      case 1
         centerCode = [0,1];
      case 2
         centerCode = [];
   end
   Best = getPollCenter(centerCode,Cache.Filter,Cache,Options.hRho);
   if length(Best) == 1
      [BestF,BestI] = deal(Best);  
   else
      [BestF,BestI] = deal(Best(1),Best(2));
   end
   fplusxi = BestF.f + RunData.fxi;
   hplusxi = min(BestI.h + RunData.hxi,Cache.Filter.hmax);
   ePollF = find([N.h] <= Options.hmin & [N.f] < fplusxi);
   ePollH = find([N.h] >  Options.hmin & [N.h] < hplusxi);

   % Determine the order in which extended poll centers get polled around
   if ~isempty(ePollF)
      [~,order] = sort([N(ePollF).f]);
      ePollF = ePollF(order);
   end
   if ~isempty(ePollH)
      [~,order] = sort([N(ePollH).h]);
      ePollH = ePollH(order);
   end
   if RunData.pcenter.h <= Cache.Filter.hmin
      ePoll = [ePollF, ePollH];
   else
      ePoll = [ePollH, ePollF];
   end

   % Search/Poll around each selected discrete neighbor using local filter
   for k = 1:length(ePoll)

      % Create MVP filter and populate it with the discrete neighbor
      pcenter = N(ePoll(k));
      Filter = createFilter(Options.hmin,hplusxi,1,0,[]);
      [~,Filter] = updateFilter(pcenter,Filter,Cache);
      if Problem.adaptOmega
         [Problem.Omega,RunData.scale] = updateOmega(Problem,Options,pcenter);
      end

      % Begin EXTENDED SEARCH step using same Search type is in the SEARCH step
      if Options.mvp1Surrogate
         success = -1;
      else
         TempRunData = RunData;
         TempOptions = Options;
         TempRunData.pcenter = pcenter;
         TempOptions.Poll.complete = Options.EPoll.complete;
         for j = 1:Options.nSearches
             TempOptions.Search(j) = Options.Search(j);
         end
         Cache        = getappdata(0,Problem.nameCache);
         startCache   = Cache.size + 1;
         [success,TempFilter] = search(Problem,TempOptions,TempRunData, ...
                                       Cache.Filter);
         Cache        = getappdata(0,Problem.nameCache);
         endCache     = Cache.size;
         Cache.Filter = TempFilter;
         setappdata(0,Problem.nameCache,Cache);

      % If EXTENDED SEARCH was unsuccessful, test versus local filter
         if success < 0
            for j = startCache:endCache
               [~,Filter] = updateFilter(Cache.iterate(j),Filter,Cache);
            end
%            pcenter = getPollCenter(Options.Poll.center,Filter,Cache);
         end
      end

      % Begin EXTENDED POLL
      while success < 0
         startCache = Cache.size + 1;
         Options.Poll = Options.EPoll;
         [ePsuccess,Filter] = poll('E',Problem,Options,RunData,pcenter,Filter);
         Cache = getappdata(0,Problem.nameCache);
         endCache = Cache.size;

         % If Extended Poll was unsuccessful, try next discrete neighbor
         if ePsuccess < 0, break; end

         % If Extended Poll step was successful, determine new incumbent
         if Options.EPoll.complete
            for j = startCache:endCache
               iterate = Cache.iterate(j);
               [unfiltered,TempFilter] = updateFilter(iterate,Cache.Filter,...
                                                      Cache);
               if unfiltered
                  success = j;
                  Cache.isSurrPoint(j) = 1;
               end
               Cache.Filter = TempFilter;
               setappdata(0,Problem.nameCache,Cache);
            end
         else
            centerCode = Options.Poll.center;
            if Options.useFilter == 2, centerCode = []; end
            [iterate,iCache] = getPollCenter(centerCode,Filter,Cache,Options.hRho);
            if Options.runStochastic, Cache.size = iCache; end
            [unfiltered,TempFilter] = updateFilter(iterate,Cache.Filter,Cache);
            if unfiltered
                success = k;
                Cache.isSurrPoint(iCache) = 1;
            end
            Cache.Filter = TempFilter;
            setappdata(0,Problem.nameCache,Cache);
         end
      end
      if success > -1 && ~Options.EPoll.complete, break; end
   end
end
return

%*******************************************************************************
% update:  Updates MADS run parameters.
% ------------------------------------------------------------------------------
% Called by: mads, processInput
% Calls:     getPollCenter, updateOmega, plotHistory, getSimplexGradient
% VARIABLES:
%  Problem             = structure defining optimization problem
%    .nameCache        =   name of the base workspace Cache variable
%    .Omega            =   structure defining linear constraints
%    .isMVP            =   flag indicating an MVP problem
%    .adaptOmega       =   flag for updating an MVP's Omega parameters
%  RunData             = structure holding MADS run statistics
%    .delta            =   current poll size
%    .nIter            =   iteration count
%    .nFunc            =   cumulative number of function evaluations
%    .grad             =   cumulative number of gradient evaluations
%    .time             =   cumulative CPU time used
%    .nFails           =   current count of consecutive Poll failures
%    .nCacheHits       =   cumulative number of Cache hits
%    .nFunc0           =   function evaluations during initialization
%    .pcenter          =   current poll center
%    .iCache           =   vector of Cache indices of all previous Poll centers
%    .bSet             =   vector used by MADS in constructing Poll directions
%    .fxi              =   f-value threshold for triggering Extended Poll
%    .hxi              =   h-value threshold for triggering Extended Poll
%    .goodD            =   index of the most recent successful direction
%    .sd               =   simplex gradient at new poll center
%    .stopRun          =   flag for stopping run immediately
%    .scale            =   scale factors for each continuous variable
%  Options             = structure holding MADS parameters
%    .computeGrad      =   flag for computing any available gradients
%    .countCache       =   flag for counting Cache points as function calls
%    .delta0           =   initial poll size
%    .deltaMax         =   maximum poll size 
%    .meshRefine       =   mesh size refinement factor
%    .meshCoarsen      =   mesh size coarsening factor
%    .accelerate       =   flag for decreasing the mesh refinement factor
%    .EPoll.fTrigger   =   f-value trigger for executing extended POLL
%    .EPoll.hTrigger   =   h-value trigger for executing extended POLL
%    .pollCenter       =   code indicating which poll center to select
%    .saveHistory      =   saves poll center to text file
%    .plotHistory      =   turns on/off plot of f-value vs. #f-evals
%    .plotColor        =   string indicating color of history plot line
%    .stophandle       =   handle for external object that terminates run
%    .runUntilFeasible =   flag for running MADS only until feasible
%    .runOneIteration  =   flag for running MADS one iteration at a time
%  success             = flag indicating iteration successful/failure
%  Cache               = database of all previously computed iterates
%    .size             =   number of iterates currently in the Cache
%    .nHits            =   number of Cache hits
%  Filter              = structure containing filter data
%    .hmax             =   maximum allowable constraint violation h-value
%    .feasible         =   indices of best feasible points in the Cache
%  meshScaleFactor     = factor for increasing or decreasing mesh size
%  BestF               = best feasible iterate found thus far
%*******************************************************************************
function [Problem,RunData] = update(Problem,Options,RunData,success)

Cache = getappdata(0,Problem.nameCache);

% Update run statistics for either initial or non-initial iterates
if success == -9
   RunData.delta       = Options.delta0;
   RunData.nIter       = 0;
   RunData.nFunc0      = (Options.countCache)*Cache.size;
   RunData.nFunc       = RunData.nFunc0;
   RunData.sdev        = 0;
   RunData.grad        = Options.computeGrad;
   RunData.nPoll       = 0;
   RunData.nFails      = 0;
   RunData.meshRefine  = Options.meshRefine;
   RunData.meshCoarsen = Options.meshCoarsen;
   RunData.hxi         = Options.EPoll.hTrigger*Cache.Filter.hmax;
   RunData.iCache      = [];
   RunData.bSet        = [];
   RunData.success     = 0;
else
   meshScaleFactor    = (success >=  1)*RunData.meshCoarsen + ...
                        (success <= -1)*RunData.meshRefine  + ...
                        (success ==  0);
   RunData.delta      = min(meshScaleFactor*RunData.delta, Options.deltaMax);
   RunData.nIter      = RunData.nIter + 1;
   RunData.nFunc      = sum(Cache.nFunc(1:Cache.size)) - RunData.nFunc0;
   RunData.grad       = RunData.grad + (Options.computeGrad && success >= 0);
   RunData.nFails     = (success == -1)*(RunData.nFails + 1);
   RunData.meshRefine = RunData.meshRefine/(1+((success == -1) && Options.accelerate));
   RunData.success    = success;
end
RunData.time       = cputime - RunData.tic;
RunData.nCacheHits = Cache.nHits;

% Add to matrix of indices of poll centers
switch Options.useFilter
   case 0
      centerCode = 0;
   case 1
      centerCode = Options.Poll.center;
   case 2
      centerCode = [];
end
[RunData.pcenter,iCache] = getPollCenter(centerCode,Cache.Filter,Cache,Options.hRho);
if isempty(RunData.iCache) || (iCache ~= RunData.iCache(end))
   RunData.iCache = [RunData.iCache, iCache];
end

% Update hmax for progressive barrier approach
if success >= 0 && Options.useFilter == 2
   Cache.Filter.hmax = max([RunData.pcenter(:).h]);
end

% Construct mesh size from poll size
switch Options.Poll.type
   case {'MADS_2n','MADSr_2','MADS_2'}
      RunData.delta_m = RunData.delta^2;
   case {'MADS_n+1'}
      RunData.delta_m = (RunData.delta/RunData.pcenter.n)^2;
   case {'OrthoMADSr_2n','OrthoMADS_2n'}
      RunData.delta_m = RunData.delta^2;
   case {'OrthoMADSr_n+1','OrthoMADS_n+1'}
      RunData.delta_m = RunData.delta^2;
   otherwise
      RunData.delta_m = RunData.delta;
end

% Compute simplex gradient if it will be used for poll ordering
if strcmp(Options.Poll.order,'SimplexGradient')
   RunData.sd = getSimplexGradient(100,Cache,iCache,RunData.delta_m,success);
end

% Update StopRun flag, as appropriate
RunData.stopRun = 0;
if ~isempty(Options.stophandle)
   RunData.stopRun = get(Options.stophandle,'UserData');
end
if Options.runUntilFeasible && ~isempty(Cache.Filter.feasible)
   RunData.stopRun = 1;
end
if Options.runOneIteration || exist(Options.abortFile,'file')
   RunData.stopRun = 1;
end
if RunData.stopRun
   drawnow
end

% For MVP, update Omega, scale factors, and extended poll thresholds
if Problem.isMVP && success ~= -1
   if Problem.adaptOmega
      [Problem.Omega,RunData.scale] = updateOmega(Problem,Options, ...
                                                  RunData.pcenter);
   end
   BestF = getPollCenter(0,Cache.Filter,Cache,Options.hRho);
   RunData.fxi = Options.EPoll.fTrigger*max(Options.Term.delta,abs(BestF.f));
end
if Options.plotHistory == 2
   plotHistory(Options.hplothandle,Options.plotColor,Cache);
end
return

%*******************************************************************************
% FUNCTIONS CALLED BY MADS SEARCH ROUTINE
%*******************************************************************************
%*******************************************************************************
% lhs:  Perform a Latin Hypercube Search for integer vectors
% ------------------------------------------------------------------------------
% Called by: search
% VARIABLES:
%  z        = matrix whose columns are LH vectors on the unit hypercube
%  [LB,UB]  = lower and upper bound vectors
%  nPoints  = number of LHS points to be generated
%  nGrid    = number of bins in each row or column (>= nPoints)
%  addNoise = flag for adding noise within each LHS bin
%  plist    = cell array of allowable categrical variable values
%  n        = dimension of sample space
%  indP     = indices of categorical variable values
%*******************************************************************************
function z = lhs(LB,UB,nPoints,nGrid,addNoise,plist)

% Latin hypercube search on [0,1]
n     = length(UB);
nx    = n - length(plist);
nGrid = max(nGrid,nPoints);
z     = zeros(n,nGrid);
for k = 1:n
   z(k,:) = (randperm(nGrid) - 1 + ~~addNoise*rand(1,nGrid))/nGrid;
end
z = z(:,1:nPoints);

% Error check and map to bound constrained region
if isempty(z)
   error('mads:lhs','No points generated by LHS');
end
z = repmat(LB,1,nPoints) + diag(UB-LB)*z;
z(nx+1:n,:) = round(z(nx+1:n,:));
return

%*******************************************************************************
% ps: Performs Particle Swarm Optimization, as part of Search.
% ------------------------------------------------------------------------------
% Called by: search
% VARIABLES:
%  xBestG           = best global solution found
%  fBestG           = objective function value of the best global solution found
%  termflag         = flag indicating if velocity termination condition is met
%  f                = objective function handle
%  B                = n x 2 matrix of lower and upper bounds
%  p                = current iterate categorical variable values
%  Options          = structure of particle swarm options
%    .maxIterations =   maximum number of iterations
%    .nParticles    =   number of particles in the swarm
%    .c             =   2 constants that weigh local versus global
%    .vTol          = velocity termination tolerance
%  PS               = structure of pre-existing particle swarm data
%    .x             =   matrix of current particle positions (1 particle/column)
%    .xBest         =   matrix of best particle positions (1 particle/ column)
%    .fBest         =   vector of best particle objective function values
%  Default          = structure of particle swarm default options
%    .maxIterations =   maximum number of iterations
%    .nParticles    =   number of particles in the swarm
%    .c             =   2 constants that weigh local versus global
%    .vTol          = velocity termination tolerance
%  LB,UB            = vectors of lower and upper bounds
%  n                = problem dimension
%  fnames           = temporary storage of the Options structure field names
%  x                = matrix of particle positions (1 particle/column)
%  v                = matrix of particle velocities (1 particle/column)
%  xBest            = matrix of best particle positions (1 particle/column)
%  fBest            = vector of best particle objective function values
%  ind              = index used to enforce bound constraints
%  fx               = vector of current particle objective function values
%  iBestG           = index of particle with best global obj function value
%  xCurrent         = temporary storage of a single particle position
%  vnorm            = vector of velocity norms for each particle
%*******************************************************************************
function [xBestG,fBestG,termflag] = ps(f,B,p,Options)

% Set options and defaults
Defaults.maxIterations = 1000;
Defaults.nParticles    = 50;
Defaults.c             = [1.3,1.3];
Defaults.vTol          = 1e-4;
Defaults.nPoints       = 1;

if nargin < 3
   Options = Defaults;
else
   fnames = fieldnames(Defaults);
   for k = 1:length(fnames)
      if ~isfield(Options, fnames{k})
         Options.(fnames{k}) = Defaults.(fnames{k});
      end
   end
end

%Retrieve pre-existing particle swarm, if available
if isappdata(0,'PS')
   PS = getappdata(0,'PS');
end

%Initialize particle positions and velocities
n  = size(B,1);
if isappdata(0,'PS') && size(PS.xBest,1) == n
   x     = PS.x;
   v     = PS.v;
   xBest = PS.xBest;
   fBest = PS.fBest;
else
   x = zeros(n,Options.nParticles);
   v = zeros(n,Options.nParticles);
   for k = 1:Options.nParticles
      x(:,k) = B(:,1) + (B(:,2) - B(:,1)).*rand(n,1);
      v(:,k) = rand(n,1);
   end
   xBest = x;
   fBest = Inf*ones(Options.nParticles,1);
end

for i = 1:Options.maxIterations
   fx = zeros(1,Options.nParticles);
   for k = 1:Options.nParticles

      % Enforce bound constraints
      ind = find(x(:,k) < B(:,1));
      x(ind,k) = B(ind,1);
      ind = find(x(:,k) > B(:,2));
      x(ind,k) = B(ind,2);

      % Evaluate objective function at a particle position
      switch nargin(f)
      case {1}
         fx(k) = feval(f,x(:,k));
      case{2}
         fx(k) = feval(f,x(:,k),p);
      otherwise
         error('PS:FuncEval','Invalid number of input arguments');
      end

      % Store best particle position and function value
      if fx(k) < fBest(k)
         fBest(k)   = fx(k);
         xBest(:,k) = x(:,k);
      end
   end

   % Get best global solution
   [fBestG,iBestG] = sort(fx);
   xBestG          = x(:,iBestG);

   % Update particle positions and velocities
   vnorm = zeros(1,Options.nParticles);
   for k = 1:Options.nParticles
      xCurrent = x(:,k);
      v(:,k)   = v(:,k) + Options.c(1)*rand*(xBest(:,k)  - xCurrent) ...
                        + Options.c(2)*rand*(xBestG(:,1) - xCurrent);
      x(:,k)   = xCurrent + v(:,k);
      vnorm(k) = norm(v);
   end

   % Termination test
   if max(vnorm) < Options.vTol
      break
   end
end

% Store PS data
PS.x     = x;
PS.v     = v;
PS.xBest = xBest;
PS.fBest = fBest;
setappdata(0,'PS',PS);

% Return the best nPoints solutions found
nPoints  = min(Options.nPoints,Options.nParticles);
xBestG   = xBestG(:,1:nPoints);
fBestG   = fBestG(1:nPoints);
termflag = max(vnorm) < Options.vTol;
return

%*******************************************************************************
% ea:  Perform Evolutionary Algorithm on a mesh, as part of Search.
% ------------------------------------------------------------------------------
% Called by: search, optimizeSurrogate
% Calls:     setPenaltyTerm
% VARIABLES:
%  X         = set of good points to be evaluated
%  fX        = function values of points in X
%  Problem   = structure of optimization problem data
%  Options   = structure of user options
%  nPoints   = number of points to return
%  pcenter   = current Poll center
%  iterate   = vector of iterates
%  Filter    = index into iterate of filter and best feasible points
%  arg       = temporary string used in writing penalty function file
%  n         = number of continuous variables
%  a         = parameter used in penalty function used in surrogate
%  f         = name of penalty function
%  fid       = file ID (handle) of penalty function
%  LB,UB     = storage of upper and lower variable bounds for EA algorithm use
%  options   = options for EA optimizer
%  es,result = GA optimizer output variables
%*******************************************************************************
function [X,fX,nFunc] = ea(f,LB,UB,Options,Search)

% Set up variable bounds for cmaes (requires full set of finite variable bounds)
ind = ~isfinite([LB,UB]);
LB(ind(:,1)) = -1/eps;
UB(ind(:,2)) =  1/eps;

% Set up initial point and sigma value
x0                = (UB + LB)./2;
x0(ind(:,1))      = UB(ind(:,1));
x0(ind(:,2))      = LB(ind(:,2));
x0(all(ind,2))    = 0;
sigma             = (UB - LB)./3;
sigma(ind(:,1))   = min(1,UB(ind(:,1)));
sigma(ind(:,2))   = min(1,LB(ind(:,2)));
sigma(all(ind,2)) = 1;

% Set up appropriate option values
Search.cmaes.Stop.maxFunEval = 50000;
Search.cmaes.Stop.tolX       = Options.Term.delta;
Search.cmaes.lBounds         = LB;
Search.cmaes.uBounds         = UB;
Search.cmaes.debug           = Options.debug;

% Call cmaes solver enough to generate nPoints minimizers
x    = x0;
X    = zeros(length(UB),Search.nPoints);
fmin = zeros(1,Search.nPoints);
for k = 1:Search.nPoints
   [xmin,fmin(k),nFunc] = cmaes_mod(f,x,sigma,Search.cmaes);
   X(:,k) = xmin;
   x      = xmin;
end
[fX, order] = sort(fmin);
X = X(:,order);
return

%*******************************************************************************
% gaPenalty:  Penalty function used by cmaes.
% ------------------------------------------------------------------------------
% Called by: ea
% VARIABLES:
%  z       = penalty function value at a point x
%  x       = point to be evaluated
%  p       = categorical variable values of point to be evaluated
%  penalty = penalty for infeasibility
%  f       = function value at point x
%  fx      = objective function value at point x
%  cx      = vector of constraint function values at point x
%  A       = linear constraint coefficient matrix
%  l       = vector of lower bounds of linear constraints
%  u       = vector of upper bounds of linear constraints
%  b       = vector of reformulated linear constraint bounds
%*******************************************************************************
function z = eaPenalty(x,p,f,A,lb,ub,penalty)
if isempty(p)
   [fx,cx] = feval(f,x);
else
   [fx,cx] = feval(f,x,p);
end
if ~isempty(A)
   n        = length(x);
   A        = [ A(n+1:end,:);  -A(n+1:end,:)];
   b        = [ub(n+1:end,:); -lb(n+1:end,:)];
   ind      = ~isfinite(b);
   b(ind)   = [];
   A(ind,:) = [];
   cx       = [cx; A*x - b];
end
z = fx + penalty*norm((cx > 0).*cx)^2;

return

%*******************************************************************************
% getDataSites:  Set up the data sites for constructing an SPS surrogate.
% ------------------------------------------------------------------------------
% Called by: search
% Calls: cat2cont
% VARIABLES:
%  Sites          = vector of iterates to be used as data sites
%    .x           =   vector of continuous variables
%    .f           =   objective function value of .x
%    .sf          =   surrogate objective function values of .x
%  Problem        = structure of optimization problem parameters
%    .isMVP       = flag indicating an MVP problem is being solved
%    .Omega.plist = cell array of lists of possible categorical variable values
%  Cache          = database of all previously computed iterates
%  Search         = structure containing Search parameters 
%    .file        = handle of the user surrogate file
%    .type        = string indicating the type of Search being performed
%  mvp1Surrogate  = flag for using 1 surrogate instead of many for MVP problems
%  debug          = numeric code for printing debugging lines
%*******************************************************************************
function Sites = getDataSites(Problem,Cache,Search,mvp1Surrogate,digits,debug)

% Get rid of data sites with infinite function values
Cache = Cache(isfinite([Cache.f]));
Sites = Cache;

% If MVP, append categorical variables to the end of the continuous variables
if Problem.isMVP && mvp1Surrogate
   for k = 1:length(Sites)
      Sites(k) = cat2cont(Sites(k),Problem.Omega.plist);
   end
end

% Evaluate the surrogate at each previously unevaluated data site
if startsWith(Search.type,'SPS') || strcmpi(Search.type,'CUSTOMS')
   for k = 1:length(Sites)
      if isempty(Sites(k).sf)
         switch nargin(Search.file)
            case {1}
               if nargout(Search.file) == 1
                  Sites(k).sf = feval(Search.file,Sites(k).x);
                  Sites(k).sc = [];
               else
                  [Sites(k).sf, Sites(k).sc] = feval(Search.file,Sites(k).x);
               end
            case {2}
               if nargout(Search.file) == 1
                  Sites(k).sf = feval(Search.file,Cache(k).x,Cache(k).p);
                  Sites(k).sc = [];
               else
                  [Sites(k).sf, Sites(k).sc] = ...
                                feval(Search.file,Cache(k).x,Cache(k).p);
               end
            otherwise
               error('mads:getDataSites',...
                     'Surrogate file must have 1 or 2 input arguments');
         end
      end
      if debug >= 3
         disp(' ');
         disp(['   True objective f-value = ', num2str(Sites(k).f, digits)]);
         disp(['   User surrogate f-value = ', num2str(Sites(k).sf,digits)]);
      end

      % For SPS surrogate, apply surrogate to difference between SPS and truth
      Sites(k).f = Sites(k).f - Sites(k).sf;
      Sites(k).c = Sites(k).c - Sites(k).sc;
      if debug >= 3
         disp(['   Additive SPS f-value   = ', num2str(Sites(k).f,digits)]);
         disp(' ');
      end
   end
end

% Check for sufficient number of data sites
nSites = length(Sites);
if ~nSites
   error('getDataSites:nSites','No data sites found');
end

return

%*******************************************************************************
% daceSurrogate:  Set up a DACE surrogate.
% ------------------------------------------------------------------------------
% Called by: search
% Calls:     constructSurrogate, updateSurrogateData, getCondR, normDL, 
%            recalKrigingSurrogate, getLowerBoundMLE, getUpperBoundMLE
% VARIABLES:
%  surrogate       = structure containing surrogate information
%    .recalibrator =   name of the surrogate recalibration function
%    .evaluator    =   name of the surrogate evaluation function
%    .merit        =   surrogate merit function parameter
%    .local        =   flag indicating if surrogate is restricted to a region
%    .param        =   flag specifying something other than objective be used
%    .searchFile   = name of user-provided surrogate file
%  Search          = structure of parameters for the current Search step
%    .dace         =   substructure of parameters for the DACE surrogate
%      .reg        =     regression function handle
%      .corr       =     corrlation function handle
%      .theta      =     correlation parameter that recalSurrogate optimizes
%      .lower      =     lower bound for theta
%      .upper      =     upper bound for theta
%      .isotropic  =     flag indicating is all theta values are the same
%  iterate         = vector of iterates (all in the Cache)
%  Filter          = indices into iterate that contain the filter
%  optimizer       = string identifying the optimizer for the surrogate
%  nIter           = current number of iterations
%  n               = number of points used to build surrogate
%*******************************************************************************
function surrogate = daceSurrogate(Search,Sites,Filter,Omega,nIter,debug)

Search.nc  = length(Sites(end).c);
Search.n   = Sites(end).n;
nResponses = Search.nc + 1;
order      = str2double(Search.dace.reg(end));
if order > 0 && nSites < prod(Search.n+(1:order))/order
   warning('dace:nSites','Not enough data sites; using regpoly0 instead.');
   Search.dace.reg = 'regpoly0';   
end

% Set the number of theta parameters (1 or 2 for isotropic, n or n+1 if not)
isExpG = strcmp(Search.dace.corr,'correxpg');
n      = isExpG + 1 + (Search.n - 1)*(1 - Search.dace.isotropic);

% Update surrogates with new sites and set recal if ill-conditioned correlation
if isappdata(0,'SUR')
   surrogate = getappdata(0,'SUR');
   for i = 1:length(surrogate)
      surrogate(i) = recalKrigingSurrogate(Sites,surrogate(i),Filter);
   end
   testKrigingInterpError('Sites',surrogate);
   condR = nan(nResponses,1);
   for k = 1:nResponses
      condR(k) = surrogate(end).(surrogate(end).resp{k}).condR;
   end
   nGP   = length(surrogate);
   recal = any(condR > Search.dace.maxCondR);
   arg   = {surrogate(end)};
else
   nGP   = 0;
   recal = 1;
   arg   = {};
end

% Start full recalibration
if recal
   if debug >= 2
      fprintf('\n   Calibrating DACE Surrogate ...\n');
   end

   % Initialize surrogate, parse data sites, get responses, set bounds on theta
   surrogate0 = constructSurrogate(Search,Omega,nIter);
   surrogate0 = updateSurrogateData(Sites,surrogate0,Filter);
   if surrogate0.usePenalty, nResponses = 1; end
   lower = nan(n,nResponses);
   upper = nan(n,nResponses);
   for k = 1:nResponses
      f  = @(x) getCondR(x,surrogate0,k) - Search.dace.maxCondR;
      g  = @(x) normDL(x,surrogate0,k, Search.dace.mleTolFinDiff) ...
                                     - Search.dace.minNormDL;
      switch nGP
         case 0
            x0 = eps;
            lo = fzero(f,x0);
            x0 = lo;
            hi = getUpperBoundMLE(g,x0,Search.dace.mleTolBounds);
         case 1
            x0 = eps;
            lo = fzero(f,x0);
            hi = 1/minDistance(Sites);
         case 2
            lo = surrogate(end).(surrogate(end).resp{k}).theta;
            hi = 1/minDistance(Sites);
      end
      lower(:,k) = lo.*ones(n,1);
      upper(:,k) = hi.*ones(n,1);
   end
   theta0     = (lower + upper)/2;
   surrogate0 = recalKrigingSurrogate(Sites,surrogate0,Filter,theta0, ...
                                      lower,upper,Search.dace.mleOpt,arg{:});

   % Add new surrogate to current array of surrogates
   iGP = 1 + (nGP > 0);
   surrogate(iGP) = surrogate0;
   testKrigingInterpError('Full',surrogate);

   % Print debugging information
   if debug >= 2
      fprintf('   Number of Gaussian processes:  %d\n', length(surrogate));
      fmt1 = '   %s-theta(%3d): Value = %25.16f, Bounds = [%25.16f,%25.16f]\n';
      fmt2 = '   %s-lambda:     Value = %25.16f, Bounds = [%25.16f,%25.16f]\n\n';
      for k = 1:nResponses
         resp  = surrogate0.resp{k};
         theta = surrogate0.(resp).theta;
         for j = 1:length(theta)
            fprintf(fmt1,resp,j,theta(j),lower(j,k),upper(j,k));
         end
         if surrogate0.(resp).lambda
            fprintf(fmt2,resp,surrogate0.(resp).lambda,0.5,1);
         end
      end
   end
end
setappdata(0,'SUR',surrogate);
return

%*******************************************************************************
% getCondR:  Approximate the condition number of the correlation matrix
% ------------------------------------------------------------------------------
% Called by: daceSurrogate
% Calls:     dacefit (DACE)
% VARIABLES:
%  condR     = approximate condition number of the correlation matrix 
%  theta     = kriging parameters that need to be fitted
%  surr      = structure of surrogate information
%    .X      =   continuous variables values of data sites
%    .Y      =   response values for each data site
%    .reg    =   function handle for kriging regression  function
%    .corr   =   function handle for kriging correlation function
%  dmodel    = kriging model structure
%    .condR  = correlation matrix condition number
%*******************************************************************************
function condR = getCondR(theta,surr,k)
theta  = max(eps,theta);
dmodel = dacefit(surr.X,surr.Y(:,k),surr.arg{:},theta);
condR  = dmodel.condR;
return

%*******************************************************************************
% normDL:  Compute the norm of the derivative of the log-likelihood function.
% ------------------------------------------------------------------------------
% Called by: daceSurrogate
% Calls:     dacefit (DACE)
% VARIABLES:
%  dLnorm         = norm of approximate derivative of log-likelihood function 
%  theta          = kriging parameters that need to be fitted
%  surr           = structure of surrogate information
%    .X           = continuous variables values of data sites
%    .Y           = response values for each data site
%    .reg         = function handle for kriging regression  function
%    .corr        = function handle for kriging correlation function
%  dmodel         = kriging stucture at h
%    .gamma       =   kriging parameter
%    .C           =   Cholesky factor of the correlation matrix
%  h              = tolerance used to approximate derivative
%  dmodel2        = kriging structure at theta+h
%  C              = Cholesky factor of correlation matrix
%  R1/R2          = correlation matrices at theta and theta+h
%  dR             = approximate derivative of correlation
%  dL             = approximate derivative of log-likelihood function
%*******************************************************************************
function dLnorm = normDL(theta,surr,k,h)
dmodel  = dacefit(surr.X,surr.Y(:,k),surr.arg{:},theta);
dmodel2 = dacefit(surr.X,surr.Y(:,k),surr.arg{:},theta+h);
Z       = dmodel.gamma';
C       = dmodel.C;
R1      = C*C';
R2      = dmodel2.C*dmodel2.C';
dR      = (R2 - R1)./h;
dL      = ((Z'*dR*Z)/dmodel.sigma2 - trace(C'\(C\dR)))/2;
dLnorm  = norm(dL);
return

%*******************************************************************************
% getUpperBoundsMLE:  Get the upper bound for an MLE optimization
% ------------------------------------------------------------------------------
% Called by: daceSurrogate
% Calls:     bisect
% VARIABLES:
%  upper = upper bound on theta
%  f     = function handle used for computing bound
%  theta = variable for which bound is computed
%  tol   = structure containing surrogate information at theta
%  fx    = temporary storage of f(theta)
%*******************************************************************************
function upper = getUpperBoundMLE(f,theta,tol)

fx = f(theta);
if fx > 0
   while fx > 0
      theta = 2*theta;
      fx    = f(theta);
    end
   upper = bisect(f,theta/2,theta,tol,1);
else
   upper = theta;
end

return

%*******************************************************************************
% minDistance:  Minimum L-inf distance between any two iterates in a set.
% ------------------------------------------------------------------------------
% Called by: daceSurrogate
% VARIABLES:
%  minD    = minimum distance
%  iterate = vector of iterates (.x/.p = continuous/categorical variable values)
%  n       = number of iterates
%  dist    = distance between two specific iterates
%*******************************************************************************
function minD = minDistance(iterate)
minD = inf;
n  = length(iterate);
%nx = length(iterate(1).x);
%np = length(iterate(1).p);
%f  = @(x,y) ~isequal(x,y);
for i = 1:n
   for j = i+1:n
%       dx = abs(iterate(i).x - iterate(j).x);
%       dx = [dx; zeros(np,1)];
%       if ~isempty(iterate(i).p)
%          dx(nx+1:end) = cellfun(f,iterate(i).p,iterate(j).p,'UniformOutput',1);
%       end
      minD = min(minD,max(abs(iterate(i).x - iterate(j).x)));
   end
end
return

%*******************************************************************************
% constructSurrogate:  Set up a surrogate structure.
% ------------------------------------------------------------------------------
% Called by: daceSurrogate, nwSurrogate, rbfSurrogate, customSurrogate
% VARIABLES:
%  surrogate       = structure containing surrogate information
%    .recalibrator =   name of the surrogate recalibration function
%    .evaluator    =   name of the surrogate evaluation function
%    .param        =   cell array of parameters needed by recalibrator
%    .isotropic    =   flag indicating an isotropic DACE correlation model
%    .optimizer    =   name of surrogate optimizer
%    .usePenalty   =   flag for using penalty function to handle constraints
%    .merit        =   surrogate merit function parameter
%    .local        =   flag indicating if surrogate is restricted to a region
%    .param        =   flag specifying something other than objective be used
%    .searchFile   = name of user-provided surrogate file
%  Search          = structure of parameters for the current Search step
%  nIter           = current number of iterations
%  recalibrator    = name of file used to recalibrate the surrogate
%  evaluator       = name of file used to evaluate the surrogate at a point
%*******************************************************************************
function surrogate = constructSurrogate(Search,Omega,nIter)

surrogate.resp = [{'f'}, arrayfun(@(k)['c',int2str(k)],1:Search.nc, ...
                                       'UniformOutput',0)];
switch upper(Search.type)
   case {'DACE','SPS-DACE'}
      surrogate.recalibrator = 'dacefit';
      surrogate.evaluator    = 'predictor';
      surrogate.isotropic    = Search.dace.isotropic;
      surrogate.arg          = {str2func(Search.dace.reg), ...
                                str2func(Search.dace.corr)};
   case {'RBF','SPS-RBF'}
      surrogate.recalibrator = 'buildRBF';
      surrogate.evaluator    = 'evalRBF';
      surrogate.param        = {Search.rbf.kernel,Search.rbf.poly};
   case {'NW','SPS-NW'}
      surrogate.recalibrator = 'buildNW';
      surrogate.evaluator    = 'evalNW';
      surrogate.arg          = {Search.nw.kernel,Search.nw.sigma, ...
                                Search.nw.lower, Search.nw.upper};
   case 'CUSTOM'
      [recalibrator,evaluator,arg] = feval(Search.file,Search.n);
      surrogate.recalibrator = recalibrator;
      surrogate.evaluator    = evaluator;
      surrogate.arg          = arg;
   otherwise
end
surrogate.optimizer    = Search.optimizer;
surrogate.usePenalty   = strcmp(surrogate.optimizer,'cmaes');
surrogate.merit        = Search.merit/(2^nIter);
surrogate.local        = Search.local;
surrogate.param        = Search.param;
surrogate.searchFile   = Search.file;
surrogate.nc           = Search.nc;
if isfield(Omega,'plist')
   surrogate.plist = Omega.plist;
end
return

%*******************************************************************************
% testKrigingInterpError:  Warn against a surrogate not fitting data
% ------------------------------------------------------------------------------
% Called by: daceSurrogate
% VARIABLES:
%   typeRecal  = string indicating type of recalibration (Full or Sites)
%   surrogate  = array of surrogate structures
%     .resp    =   cell array of response names
%     .X       =   matrix of data sites used to build kriging model
%     .Y       =   matrix of response values at each data site
%   surrVal    = array of surrogate-predicted values at each data site
%   maxErr     = maximum absolute error between truth and surrogate at data site
%   msg1/msg2  = strings used in warning message
%*******************************************************************************
function testKrigingInterpError(typeRecal,surrogate)

fprintf('   Maximum Surrogate Error (%s):\n', typeRecal);
for i = 1:length(surrogate)
   resp = surrogate(i).resp;
   for k = 1:length(resp)
      surrVal = ones(size(surrogate(i).X,1),1);
      for j = 1:size(surrogate(i).X,1)
         surrVal(j) = predictor(surrogate(i).X(j,:),surrogate(i).(resp{k}));
      end
      maxErr = max(abs(surrVal - surrogate(i).Y(:,k)));
      fprintf('      Surrogate %d (%2s)  = %23.15f\n',i,resp{k},maxErr);
      if maxErr > 1e-4
         warning on
         msg1 = ['Kriging interpolation error for response ',resp{k}];
         msg2 = [' exceeds tolerance (',num2str(maxErr),').'];
         warning('mads:recalKrigingSurrogate',[msg1,msg2]);
         warning off
      end
   end
end
return

%*******************************************************************************
% recalKrigingSurrogate:  Recalibrate a kriging surrogate.
% ------------------------------------------------------------------------------
% Called by: daceSurrogate
% Calls:     updateSurrogateData, dacefit (DACE)
% VARIABLES:
%  Site         = vector of data sites (iterates in the Cache)
%  surrogate    = structure containing surrogate information
%    .resp      =   substructure of response labels
%    .optimizer =   string identifying the optimizer used to find S
%  Filter       = indices into iterate that contain the filter
%  usePenalty   = flag for using a penalty function in surrogate optimization
%  X            = array of continuous variables values of data sites
%  Y            = array of response values of data sites
%  F,C          = objective and constraint values of data sites
%  resp         = temporary storage of response labels
%  nResponses   = number of responses
%  theta        = kriging parameters to be fitted during recalibration
%  lower/upper  = lower and upper bounds on theta for MLE fit
%  mleOpt       = structure of MLE optimization parameters
%  surr         = copy of all previous surrogate information
%  arg          = cell array of input arguments for dacefit
%*******************************************************************************
function surrogate = recalKrigingSurrogate(Sites,surrogate,Filter,varargin)

% Update surrogate data
surrogate  = updateSurrogateData(Sites,surrogate,Filter);
resp       = surrogate.resp;

% Process optional input arguments
switch nargin
   case 4
      theta = varargin{1};
   case 7
      [theta,lower,upper,mleOpt]      = deal(varargin{:});
   case 8
      [theta,lower,upper,mleOpt,surr] = deal(varargin{:});
end

% Recalibrate the surrogate for each response
for k = 1:length(resp)
   switch nargin
      case 3
         arg = {surrogate.(resp{k}).theta};
      case 4
         arg = {theta(:,k)};
      case 7
         arg = {theta(:,k),lower(:,k),upper(:,k),mleOpt};
      case 8
         arg = {theta(:,k),lower(:,k),upper(:,k),mleOpt,[surr.(resp{k})]};
      otherwise
         error('mads:recalSurrogate','Invalid number of input arguments');
   end
   surrogate.(resp{k}) = dacefit(surrogate.X,surrogate.Y(:,k), ...
                                 surrogate.arg{:},arg{:});
   if surrogate.isotropic && length(surrogate.(resp{k}).theta) > 1
      surrogate.(resp{k}).theta(2:end) = [];
   end
end

return

%*******************************************************************************
% updateSurrogateData:  Get data matrices used in kriging model construction.
% ------------------------------------------------------------------------------
% Called by: daceSurrogate,recalSurrogate,recalKrigingSurrogate
% Calls:     setPenaltyTerm
% VARIABLES:
%  Site          = vector of data sites (iterates in the Cache)
%  surrogate     = structure of surrogate information
%    .usePenalty = flag indicating use of penalty function to handle constraints
%    .param      = flag indicating surrogate construction on user parameter
%  Filter        = indices into iterate that contain the filter
%  X             = continuous variables values of data sites
%  Y             = response values for each data site
%  F,C           = objective and constraint values of iterates used for surrogate
%  H             = constraint violation function values of surrogate iterates
%  nSites        = number of points used to build surrogate
%  nConstraints  = number of nonlinear constraints
%  Param         = user defined problem parameters
%*******************************************************************************
function surrogate = updateSurrogateData(Site,surrogate,Filter)

% Delete data sites without finite objective function values
if surrogate.usePenalty
   a = setPenaltyTerm(Site,Filter);
end
Site(~isfinite([Site.f])) = [];

% Set up data sites
nSites       = length(Site);
nConstraints = length(Site(end).c);
X = [Site.x]';
F = [Site.f]';
H = [Site.h]';
C = zeros(nSites,nConstraints);
if nConstraints
   for k = 1:nSites
      C(k,:) = Site(k).c;
   end
end

% Allow surrogate construction on user-specified parameter instead of objective
if isappdata(0,'PARAM')
   Param = getappdata(0,'PARAM');
   if isfield(Param,'param') && surrogate.param
      F = zeros(nSites,length(Site(1).param));
      for k = 1:nSites
         F(k,:) = Site(k).param;
      end
   end
end

% Group response data together
if surrogate.usePenalty
   Y = F + a*H;
else
   Y = [F, C];
end

% Compute max distance between all data sites
if surrogate.local
   D = zeros(nSites,nSites);
   for i = 1:nSites
      for j = 1:i-1
         D(i,j) = norm(X(i,:) - X(j,:));
      end
   end
   maxDistance = max(D(:));
else
   maxDistance = inf;
end

% Compute trust region radius and distance and merit function penalty term
surrogate.trust = maxDistance/2;
surrogate.dist  = 10*surrogate.merit*(max(F) - min(F));

% Set output
surrogate.X = X;
surrogate.Y = Y;
surrogate.F = F;
surrogate.C = C;

return

%*******************************************************************************
% setPenaltyTerm: Get penalty parameter for applying EA to constrained problem.
% ------------------------------------------------------------------------------
% Called by: search (ea), updateSurrogateData
% VARIABLES:
%  a           = parameter used in penalty function used in surrogate or EA
%  iterate     = vector of iterates
%    .f        =   objective function value
%    .h        =   constraint violation function value
%  Filter      = index into iterate of filter and best feasible points
%    .F        =   indices of filter points
%    .feasible =   indices of best feasible solutions
%  FF          = objective function values of filter points
%  HH          = constraint violation function values of filter points
%*******************************************************************************
function a = setPenaltyTerm(iterate,Filter)

% Compute penalty parameter
FF = [iterate(Filter.F).f]';
HH = [iterate(Filter.F).h]';
if isempty(Filter.feasible)
   a = 1000;
else
   if Filter.feasible(1) < length([iterate.f])
      FF = [iterate(Filter.feasible(1)).f; FF]';
      HH = [iterate(Filter.feasible(1)).h; HH]';
   end
   a  = max((FF(1:end-1)-FF(2:end))./(HH(2:end)-HH(1:end-1)));
   if isempty(a), a = 0.5; end
end
return

%*******************************************************************************
% optimizeSurrogate:  Optimize a surrogate problem.
% ------------------------------------------------------------------------------
% Called by: search
% Calls: dacefit (DACE), buildNW (NW), buildRBF (RBF)
% VARIABLES:
%  S           = set of points generated to be evaluated
%  Problem     = structure containing optimation problem data
%    .Omega    =   structure defining feasible region
%      .A      =     matrix of linear coefficients
%      .l      =     vector of lower bounds
%      .u      =     vector of upper bounds
%  pcenter     = current Poll center
%  surrogate   = structure containing surrogate information
%  optimizer   = string identifying the optimizer used to find S
%  isSPS       = flag indicating a simplified physics surrogate
%  nPoints     = number of points to be returned by the optimizer
%  tol         = desired accuracy of optimizer
%  p           = cell array containing surrogate Functions file source code
%  f           = name of the surrogate objective function file
%  C           = name of the surrogate constraints function file
%  fid         = file ID (handle) for a surrogate problem file
%  Sur         = structure of input data for running MADS on the surrogate
%  SurRunData  = RunData for MADS run on the surrogate problem
%  SurCache    = Cache for the MADS run on the surrogate problem
%  temp1,2     = temporary storage
%  indFeasible = indices of feasible points to be evaluated
%  indFilter   = indices of infeasible points to be evaluated
%  n           = number of continuous variables
%  es_options  = options for EA optimizer
%  es,result   = EA optimizer output variables
%  sf          = surrogate function values of points from EA optimizer
%  order       = indices of EA-produced iterates, sorted by sf
%  y           = temporary storage
%*******************************************************************************
function [S,termflag] = optimizeSurrogate(Problem,pcenter,surrogate,Options, ...
                                          Search,tol)

termflag = 0;
% if isfield(Problem.Omega,'plist')
%    surrogate.plist = Problem.Omega.plist;
% end
% setappdata(0,'SUR',surrogate);

% Merge continuous and categorical variables if only one MVP surrogate
if Problem.isMVP && Options.mvp1Surrogate && ~strcmp(surrogate(end).optimizer,'mads')
   np = length(pcenter.p);
   iterate = cat2cont(pcenter,Problem.Omega.plist);
else
   np = 0;
   iterate = pcenter;
end
nc = length(pcenter.c);

% Optimize based on user-specified optimization scheme
switch surrogate(end).optimizer

% Call a Matlab optimizer
case {'fmincon','MultiStart','GlobalSearch'}

   % Make sure the Optimization Toolbox is available
   if ~license('test','Optimization_Toolbox')
       error('mads:missingToolbox', ...
             'Optimization Toolbox required (optimizeSurrogate).');
   end

   % Set options and surrogate function names
   Options = optimoptions('fmincon','MaxIterations',500, ...
                    'MaxFunctionEvaluations',5000,'Display','off','TolX',tol);
   f = @SurObj;
   C = [];
   if nc
      C = @SurCon;
   end

   % Call MultiStart or fmincon optimizer
   warning('off','all');
   if isempty(Problem.Omega.A) && isempty(C)
      Prob = createOptimProblem('fmincon','objective',f,'x0',iterate.x, ...
                                'options',Options);
      if ~strcmp(surrogate(end).optimizer,'fmincon')
         error('mads:optimizeSurrogate','MultiStart must have variable bounds');
      end
      S.x = fminunc(Prob);
   else
      n   = pcenter.n;
      nL  = size(Problem.Omega.A,1) - n;
      A   = [Problem.Omega.A(n+1:end,:), zeros(nL,np); zeros(np,n), eye(np)];       
      b   = [Problem.Omega.u(n+1:end);  iterate.x(end-np+1:end); ...
            -Problem.Omega.l(n+1:end); -iterate.x(end-np+1:end) ];
      LB   = Problem.Omega.l(1:n);
      UB   = Problem.Omega.u(1:n);
      Prob = createOptimProblem('fmincon','objective',f,'x0',iterate.x,...
                                'Aineq',[A;-A],'bineq',b,'lb',LB,'ub',UB, ...
                                'nonlcon',C,'options',Options);
      rng default;
      switch surrogate(end).optimizer
         case 'fmincon'
            S.x = fmincon(Prob);
         case 'MultiStart'
            ms  = MultiStart;
            S.x = run(ms,Prob,20);
         case 'GlobalSearch'
            gs  = MultiStart;
            S.x = run(gs,Prob);
      end
      if Problem.isMVP && Options.mvp1Surrogate
         S.x = S.x(1:length(pcenter.x));
      end
   end
   warning('on','all');
   S.p = pcenter.p;

% Call optimizer MADS
case 'mads'

   % Setup options for Surrogate optimization by MADS
   Sur.Defaults               = mads_defaults('SURROGATE');
   Sur.Problem                = Problem;
   Sur.Problem.nameCache      = Sur.Defaults.nameCache;
   Sur.Problem.typeProblem    = Sur.Defaults.typeProblem;
   Sur.Problem.File.F         = 'SurFun';
   Sur.Problem.File.C         = [Problem.File.C, '_Sur_Cache'];
   Sur.Problem.File.X         = Problem.File.X;
   Sur.Problem.fType          = 'M';
   Sur.Options                = Options.Sur;
   Sur.Options.runStochastic  = Options.runStochastic;
   Sur.Options.fixCategorical = Options.fixCategorical;
   Sur.Options.delta0         = Options.delta0;
   Sur.Options.deltaMax       = Options.deltaMax;

   % If MVP and 1 surrogate, write its Omega file if it doesn't already exist
   if Problem.isMVP && Options.mvp1Surrogate
      SurFields = {'sf','sc','sh','sgradf','sgradc','sgradh'};
      ind = isfield(iterate,SurFields);
      Sur.Problem.iterate0 = rmfield(iterate,SurFields(ind));
      Sur.Problem.isMVP    = 1;
      Sur.Options.scale    = 0;
      Sur.Options.removeRedundancy = 0;
   end

   % Pass function to MADS optimizer
   [~,~,sStats,SurCache] = mads(Sur.Problem,iterate,Sur.Options);
   termflag = sStats.stopRun;

   % Retrieve and evaluate a number of points
%    indFeasible = 1:min(Search.nPoints,length(SurCache.Filter.feasible));
%    indFilter   = 1:min(Search.nPoints,length(SurCache.Filter.F));
%    surrPoints  = [SurCache.iterate(SurCache.Filter.feasible(indFeasible)), ...
%                   SurCache.iterate(SurCache.Filter.F(indFilter))];
   S1    = SurCache.iterate(SurCache.Filter.feasible);
   S2    = SurCache.iterate(SurCache.Filter.F);
   SS    = [screenClosedConstraints(S1,Problem.File.X,Search.nPoints);
            screenClosedConstraints(S2,Problem.File.X,Search.nPoints)];
   [S.x] = SS.x;
   [S.p] = SS.p;

% Call particle swarm function as the surrogate optimizer
case 'ps'
   Options.removeRedundancy = 0;
   Options.scale = 0;
   [Problem.Omega] = updateOmega(Problem,Options,pcenter);
   LB = [ Problem.Omega.l(1:pcenter.n); iterate.x(end-np+1:end)];
   UB = [ Problem.Omega.u(1:pcenter.n); iterate.x(end-np+1:end)];
   LB(~isfinite(LB)) = -realmax;
   UB(~isfinite(UB)) =  realmax;
   Options.nPoints = 10*Search.nPoints;
   [x,~,termflag] = ps(Problem.File.F,[LB,UB],pcenter.p,Options);
   S = dealIterate(pcenter,Options.nPoints,x);
   S = screenClosedConstraints(S,Problem.File.X,nPoints);

% Call evolutionary algorithm package as the surrogate optimizer
case 'cmaes'
  Options.removeRedundancy = 0;
  Options.scale = 0;
  Problem.Omega = updateOmega(Problem,Options,pcenter);
  LB = [ Problem.Omega.l(1:pcenter.n); iterate.x(end-np+1:end)-tol ];
  UB = [ Problem.Omega.u(1:pcenter.n); iterate.x(end-np+1:end)+tol ];
  f  = @(x) feval(surrogate(end).evaluator,x,surrogate(end).f);
  X  = ea(f,LB,UB,Options,Search);
%    ind = ~isfinite([LB,UB]);
%    LB(ind(:,1)) = -1/eps;
%    UB(ind(:,2)) =  1/eps;
% 
%    % Set up initial point and sigma value
%    x0                = (UB + LB)./2;
%    x0(ind(:,1))      = UB(ind(:,1));
%    x0(ind(:,2))      = LB(ind(:,2));
%    x0(all(ind,2))    = 0;
%    sigma             = (UB - LB)./3;
%    sigma(ind(:,1))   = min(1,UB(ind(:,1)));
%    sigma(ind(:,2))   = min(1,LB(ind(:,2)));
%    sigma(all(ind,2)) = 1;
% 
%    % Set up appropriate option values
%    Search.cmaes.maxFunEval = 50000;
%    Search.cmaes.tolX       = Options.Term.delta;
%    Search.cmaes.lBounds    = LB;
%    Search.cmaes.uBounds    = UB;
%    Search.cmaes.debug      = Options.debug;
% 
%    % Call cmaes solver enough to generate nPoints minimizers
%    x  = x0;
%    X  = zeros(iterate.n,Search.nPoints);
%    sf = zeros(1,Search.nPoints);
%    for k = 1:Search.nPoints
%       [xmin,sf(k)] = cmaes_mod(surrogate(end).evaluator,x,sigma, ...
%                                Search.cmaes,surrogate(end).f);
%       X(:,k) = xmin;
%       x = xmin;
%    end
%    [~,order] = sort(sf);
%    X = X(1:pcenter.n,order);
  X = X(1:pcenter.n,:);
  S = dealIterate(pcenter,size(X,2),X);
  S = screenClosedConstraints(S,Problem.File.X,Search.nPoints);
end
return

%*******************************************************************************
% cblgs:  Constrained Balanced Local-Global Search
%      Algorithm based on AIAA 2000-4891 (Audet, Dennis, Moore, Booker, Frank)
% ------------------------------------------------------------------------------
% Called by: optimizeSurrogate
% Calls:     lhs, dealIterate, ExpectedImprovement
% VARIABLES:
%  S           = vector of iterates to be evaluated
%  Problem     = structure of Problem parameters
%    .Omega    =   substructure of bound/linear constraint informations
%      .l/.u   =     lower and upper bounds
%  pcenter     = current poll center
%    .f        =   best objective function values
%    .n        =   number of variables
%  surrogate   = structure of current DACE kriging model parameters/data
%    .resp     =   cell array of response names
%  Opt         = structure of options
%    .nGoal    =   number of points to be evaluated
%    .nCloud   =   number of points in initial dense cloud of candidate points
%    .maxCloud =   maximum number of cloud points
%    .tolFeas  =   tolerance on expected feasibility
%  surr        = optional argument containing full surrogate structure
%  y           = function handle for evaluator with fixed model data
%  resp        = cell array of response names
%  fmin        = best objective function value found thus far
%  cmax        = nonlinear constraint upper bound
%  n           = dimension of sample space
%  addNoise    = flag for adding noise within each LHS bin
%  LB, UB      = lower and uppeer bounds on each variable
%  nPoints     = number of LHS points to be generated
%  nGrid       = number of bins in each row or column (>= nPoints)
%  z           = matrix of unscaled LHS points
%*******************************************************************************
function S = cblgs(LB,UB,plist,pcenter,surrogate,Opt,xFile)

% Initialize
resp = surrogate.resp;
if ~isempty(plist)
   pUB = cellfun(@length,plist);
   LB  = [LB; ones(length(plist),1)];
   UB  = [UB; pUB(:)];
end

% CBLGS
for i = 1:Opt.maxIter
   nPoints = Opt.nCloud + (Opt.maxCloud - Opt.nCloud)*i/Opt.maxIter;

   % Generate LHS points and scale them within variable bounds
   X = lhs(LB,UB,nPoints,nPoints,0,plist)';

   % Get sufficiently feasible points
   EV = zeros(size(X,1),surrogate.nc);
   for k = 2:length(resp)
      y       = @(x) predictor(x,surrogate.(resp{k}),surrogate.(resp{k}));
      EV(:,k) = expectedImprovement(X,y,sqrt(surrogate.(resp{k}).sigma2),0,-1);
   end
   exph  = vecnorm(EV,inf,2);
   ind   = exph < Opt.tolFeas;
   nFeas = sum(ind);
   if nFeas >= Opt.nGoal, break, end
end

% Compute and sort expected improvement
X     = X(ind,:);
y     = @(x) predictor(x,surrogate.f,surrogate.f);
EI    = expectedImprovement(X,y,sqrt(surrogate.f.sigma2),pcenter.f,1);
nGoal = min(nFeas,Opt.nGoal);

% Grab the requested number of best EI points
if exist(xFile,'file')
   [~,ind] = sort(EI,'descend');
else
   [~,ind] = maxk(EI,nGoal);
end
S = dealIterate(pcenter,length(ind),X(ind,:)',plist);
S = screenClosedConstraints(S,xFile,nGoal);
return

%*******************************************************************************
% screenClosedContraints:  Screen an ordered set of iterates for closed
%                          constraint feasibility.
% ------------------------------------------------------------------------------
% Called by: search, cblgs, optimizeSurrogate
% VARIABLES:
%   S        = vector of iterates to be screened
%   xFile    = name of closed constraint file
%   nGoal    = number of (closed constraint) feasible points desired
%   nS       = length of S
%   feasible = vector of flags indicating feasible or infeasible
%*******************************************************************************
function S = screenClosedConstraints(S,xFile,nGoal)
if isempty(S), return, end
nS = length(S);
if exist(xFile,'file')
   feasible = false(nS,1);
   for k = 1:nS
      [~,feasible(k)] = evalClosedConstraints(xFile,S(k),0);
      if sum(feasible) >= nGoal
         break;
      end
   end
   S = S(feasible);
else
   S = S(1:min(nS,nGoal));
end
return

%*******************************************************************************
% ExpectedImprovement:  Compute expected improvement at a point
%      Algorithm based on AIAA 2000-4891 (Audet, Dennis, Moore, Booker, Frank)
% ------------------------------------------------------------------------------
% Called by: cblgs
% VARIABLES:
%  ei       = expected improvement value
%  x        = point being evaluated
%  y        = prediction function handle for fixed surrogate data
%  sigma    = standard deviation of responses
%  bound    = minimum function value found thus far or constraint bound
%  signFlag = flag containing either 1 or -1
%  normPDF  = function handle for the normal probability density function
%  cumPDF   = function handle for the cumuluative normal distribution
%  dist     = function handle for distance between predictor and target value
%*******************************************************************************
function ei = expectedImprovement(x,y,sigma,bound,signFlag)

if abs(signFlag) ~= 1
   error('mads:dace:ExpectedImprovement','Argument SIGNFLAG must be -1 or 1.');
end
normPDF = @(x) exp(-x.^2/2)/(sigma.*sqrt(2*pi));
cumPDF  = @(x) (1 + erf(x/sqrt(2)))/2;
dist    = @(x) signFlag*(bound - y(x));

ei = 0;
if sigma
   ei = dist(x).*cumPDF(dist(x)/sigma) + sigma.*normPDF(dist(x)/sigma);
end
return

%*******************************************************************************
% FUNCTIONS CALLED BY MADS POLL ROUTINE
%*******************************************************************************
%*******************************************************************************
% standDirections: Returns directions that positively span the tangent cone
%    at the current iterate, with respect to bound and linear constraints.
% ------------------------------------------------------------------------------
% Called by: poll
% Calls:     activeConstraints, getTangentCone, scaleDirections
% VARIABLES:
%  D        = matrix whose columns form the positive spanning set
%  x        = continuous variable values of the current iterate
%  Omega    = structure whose fields define Omega = {x: l <= Ax <= u}
%    .A     =   matrix of coefficients for the linear constraints
%    .l     =   vector of lower bounds
%    .u     =   vector of upper bounds
%  strategy = code indicating chosen poll strategy
%  tol      = tolerance for considering constraints active
%  scale    = scale factors for scaling the appropriate Poll directions
%  scheme   = scheme for handling degenerate linear constraints
%  success  = flag indicating previous iteration was successful
%  Ax       = A*x
%  n        = number of columns of A
%  I        = identity matrix
%  W        = underlying basis of Poll directions
%  Param    = structure of user-provided problem parameter data
%  B,N      = matrices whose columns together span the tangent cone
%*******************************************************************************
function D = standDirections(x,Omega,strategy,tol,scale,scheme,success)

% Initialize variables
Ax = Omega.A*x;
n  = length(x);

% Get user-provided Poll basis, if provided
I = eye(n);
W = I;
if isappdata(0,'PARAM') && strncmp(strategy,'Custom',6)
   Param = getappdata(0,'PARAM');
   if isfield(Param,'PollBasis'), W = Param.pollBasis; end
end
if ~((rank(W) == n) && isequal(size(W),[n,n]))
   error('mads:userPollBasis:dim', ...
         ['User Poll basis has incompatible size or is not a basis ',...
          '(standDirections).']);
end

% Compute B and N, and scale each column in N
if isequal(Omega.A,I)   % Bound constraints
   lactive = activeConstraints( Ax, Omega.l, Omega.u,tol);
   uactive = activeConstraints(-Ax,-Omega.u,-Omega.l,tol);
   active  = find(lactive | uactive);
   B = I(:,active);
   N = I(:,setdiff(1:n,active));
   if isempty(active), N = W; end
   B = scaleDirections(B,scale);
   N = scaleDirections(N,scale);
else                      % General linear constraints
   [B,N] = getTangentCone(x,Omega,tol,scale,scheme,success);
end

% Form directions that positively span the tangent cone at x
switch strategy
case {'Standard_n+1','Gradient_n+1','Custom_n+1'}
   D = [-sum(N,2), N, B, -B];
case {'Standard_2n', 'Gradient_2n', 'Custom_2n'}
   D = [N, -N, B, -B];
otherwise
   error('mads:poll:choice','Invalid Poll strategy type (standDirections).');
end
D(:,~any(D)) = [];
return

%*******************************************************************************
% madsDirections:  Returns a random set of directions that positively span
%    the tangent cone, with respect to bound and linear constraints.
% ------------------------------------------------------------------------------
% Called by: poll
% Calls:     scaleDirections
% VARIABLES:
%  D        = matrix whose columns form the positive spanning set
%  bSet     = matrix of previously used random MADS vectors
%  x        = continuous variable values of the current iterate
%  strategy = code indicating chosen poll strategy
%  delta0   = initial poll size parameter
%  delta    = poll size parameter
%  scale    = scale factors for scaling the appropriate Poll directions
%  n        = number of continuous variables
%  k        = index into bSet corresponding to the mesh size
%  m        = temporary storage of 2^k
%  b        = MADS direction that must remain the same for a specified mesh size
%  ihat     = randomly selected index used to construct b
%  ind      = temporary storage of randperm indices
%  L        = matrices used in forming the MADS directions
%  B        = basis from which to construct MADS directions
%*******************************************************************************
function [D,bSet] = madsDirections(x,strategy,delta0,delta,scale,bSet)

% Get index into set of previously used directions, based on poll size
n = length(x);
if n == 1
   strategy = 'MADS_2';
end
if strncmp(strategy,'Ortho',5)
   k = round(-log(delta/delta0)/log(8));
else
   k = round(-log(delta/delta0)/log(4));
end
m = 2^k;

% Get direction b that must remain the same for a specified step size
if k + 1 <= size(bSet,2)
   b = bSet(:,k+1);
   [~,ihat] = max(abs(b));
else
   switch strategy

   % Get b from Adjusted Halton sequence
   case {'OrthoMADS_2n','OrthoMADS_n+1','MADS_2'}
      b = adjHalton(n,n+1,-log(delta)/log(2));

   % Get random b based on random numbers on [-m+1,m-1], index, sign
   otherwise
      b = -m + ceil((2*m-1).*rand(n,1));
      ihat = ceil(n.*rand);
      b(ihat) = sign(ceil(2.*rand)-1.5)*m;
   end
   bSet(:,end+1) = b;   
end

% Construct basis matrix B matrix, according to the choice of MADS directions
switch strategy

% Construct B from b and (n-1)x(n-1) lower triangular matrix L
case {'MADS_2n','MADS_n+1'}
   L = zeros(n-1);
   for i = 1:n-1
      L(i,1:i-1) = -m + ceil((2*m-1).*rand(1,i-1));
      L(i,i)     = sign(ceil(2.*rand)-1.5)*m;
   end
   ind = randperm(n);
   ind(ind == ihat) = [];
   B(ind,:)         = L;
   B(ihat,:)        = 0;
   B(:,end+1)       = b;
   B = B(:,randperm(n));

% Construct othogonal B from b using Householder-like transformation
case {'OrthoMADSr_2n','OrthoMADSr_n+1','OrthoMADS_2n','OrthoMADS_n+1'}
   B = b'*b*eye(n) - 2*(b*b');
   
% Construct two-point poll from b
case {'MADSr_2','MADS_2'}
   B = b;
end

% Construct positive basis
switch strategy
case {'MADS_2n','OrthoMADSr_2n','OrthoMADS_2n','MADSr_2','MADS_2'}
   D = [B, -B];
case {'MADS_n+1','OrthoMADSr_n+1','OrthoMADS_n+1'}
   D = [B, -sum(B,2)];
otherwise
   error('mads:poll:choice','Invalid Poll strategy type (madsDirections).');
end
D = scaleDirections(D,scale);
return

%*******************************************************************************
% adjHalton:  Compute an adjusted Halton direction, based on Halton sequences.
% ------------------------------------------------------------------------------
% Called by: madsDirections
% VARIABLES:
%  b = adjusted Halton direction vector
%  n = length of the Halton direction vector
%  t = Halton sequence index
%  L = parameter used by OrthoMADS to control mesh size
%  p = vector of the first n prime numbers
%  s = vector representation of dec2base conversion
%  u = Halton direction vector
%*******************************************************************************
function b = adjHalton(n,t,L)

q = @(a,x,y) norm(round(a*x)) - y;

% Construct Halton direction vector u
t0 = n;
p  = primes(max(13,10*1.3*n*log(n)));
p  = p(t0+(1:n));
u  = zeros(n,1);
for k = 1:n
   if p(k) <= 36
      a = fliplr(dec2base(t,p(k),n)) - '0';
   else
      a = fliplr(dec2bigbase(t,p(k),n));
   end
   u(k) = sum(a./(p(k).^(1:length(a))));
end

% Construct adjusted Halton direction vector b
bound = 2^(abs(L)/2);
z     = 2*u - ones(n,1);
z     = z/norm(z);
lb    =  bound/sqrt(n) - 0.5;
ub    = (bound/sqrt(n) + 0.5)/min(abs(z));
alpha = bisect(@(alpha) q(alpha,z,bound),lb,ub,1e-8,-1);
b     = round(alpha*z);
return

%*******************************************************************************
% gradDirections:  Compute Poll directions for a gradient-pruned Poll.
% ------------------------------------------------------------------------------
% Called by: poll
% Calls:     activeConstraints, getTangentCone, dVector, scaleDirections, runLP
% VARIABLES:
%  D                 = set of directions generated by gradient-pruned polling
%  iterate           = current iterate
%  Omega             = structure describing feasible region
%  strategy          = gradient poll type
%  Filter            = structure of Cache indices of iterates in the filter
%    .hmax           = max allowable constraint violation of any filter point
%  delta             = poll size parameter
%  tol               = tolerence for determining if a constraint is active
%  scale             = direction scale factors
%  scheme            = scheme for handling degenerate linear constraints
%  success           = flag indicating previous iteration was successful
%  p                 = identifier for type of gradient-based direction
%  infGrad           = flags for identifying unavailable partial derivatives
%  g                 = gradient vector
%  d                 = gradient-based descent vector
%  B                 = indices of numerically binding nonlinear constraints
%  m                 = number of binding nonlinear constraints, plus one
%  tmax,LB,UB        = LP variables used for finding feasible descent direction
%  f,A,Aeq,b,beq     = LP variables used for finding feasible descent direction
%  X,FX              = solution of LP solve
%  y                 = feasible descent direction (y = part of X)
%  exitflag          = error flag for bad LP problem
%  goodMeshPoint = flag indicating a feasible mesh point has been found
%  z                 = closest mesh point to pcenter + y
%  I                 = identity matrix
%  Ax                = product of Omega.A and iterate.x
%  feasible          = flag for indicating if poll center is feasible
%  lactive           = flags for numerically active lower bounds of Ax
%  uactive           = flags for numerically active upper bounds of Ax
%  active            = flags for numerically active bound/linear constraints
%  [B,N]             = tangent cone generating directions
%  Bgood,Bbad        = partitioning of B to deal with unavailable partials
%  Ngood,Nbad        = partitioning of N to deal with unavailable partials
%  gradcB            = least squares approx of NLP tangent cone generators
%  DN                = descent vectors in N
%  NN                = [N, -N]
%  descent           = descent vectors in I
%  dhat              = sum of the chosen ascent and decent directions
%*******************************************************************************
function D = gradDirections(iterate,Omega,strategy,Filter,delta,tol,scale,...
                            scheme,success)

% Determine type of pruning vector
switch strategy
case 'Gradient_3n_L1'
   p = 1;
case 'Gradient_3n_L2'
   p = 2;
case 'Gradient_3n_LInf'
   p = Inf;
case 'Gradient_3n2n'
   p = 2;
otherwise
   error('mads:poll:choice','Invalid Poll strategy type (gradDirections).')
end

% Compute the descent directions that will be used for pruning
infGrad = ~isfinite(iterate.gradf);
g = iterate.gradf;
g(infGrad) = 0;
d = -dVector(g.*scale,p);
descent = diag(-dVector(g,Inf));

% Determine the index of binding and violated constraints
O  = Omega;
binding = find(iterate.c > -tol);
gradcB  = iterate.gradc(:,binding);
nB = length(binding);
if ~isempty(binding)

   %  If hmax is exceeded, then use d = -gradh
   if iterate.h > Filter.hmax
      infGrad    = [~isfinite(iterate.gradh), infGrad];
      g          = iterate.gradh;
      g(infGrad) = 0;
      d          = -dVector(g.*scale,p);
      descent    = diag(-dVector(g,Inf));

   % If hmax is not exceeded, solve LP to get descent direction in f and cB
   else
      t  = delta^2;
      c  = 0.05*ones(nB,1);
      AA = [1, iterate.gradf'; c, gradcB'];
      b  = zeros(nB+1,1);
      LB = [0; -Inf(iterate.n,1)];
      UB = [t;  Inf(iterate.n,1)];
      f  = [-1; zeros(iterate.n,1)];
      [X,~,exitflag] = runLP(f,AA,b,[],[],LB,UB);
      if exitflag == -9
         error('mads:missingToolbox', ...
               'Optimization Toolbox required (gradDirections)');
      end
      if exitflag <= 0, D = []; return, end
      y = X(2:iterate.n+1);

      % Modify descent direction so that the resulting point is on the mesh
      m = 1;
      goodMeshPoint = 0;
      while ~goodMeshPoint
         m = m + 1;
         z = delta*round((m*y - iterate.x)./(scale*delta));
         goodMeshPoint = iterate.gradf'*z < 0 && all(gradcB'*z < 0);
         if (m > 1e+3)
            d(:,end+1) = z; %#ok
            break;
         end
      end
      z = z/norm(z,p);
      if goodMeshPoint, d(:,end+1) = z; end

      % Add linearized constraints for tangent cone generation
      O.A = [O.A;  gradcB'];
      O.l = [O.l; -inf(nB,1)];
      O.u = [O.u; -iterate.c(:)];
   end
end

% Add tangent cone generators and compensate for unavailable derivatives
I = eye(iterate.n);
if isequal(O.A,I)
   Ax = O.A*iterate.x;
   lactive = activeConstraints( Ax, O.l, O.u,tol);
   uactive = activeConstraints(-Ax,-O.u,-O.l,tol);
   active  = lactive | uactive;
   Bgood   = (~infGrad(:,end) & xor(lactive,uactive));
   Ngood   = (~infGrad(:,end) & ~active);
   Bbad    = ( infGrad(:,end) &  active);
   Nbad    = ( infGrad(:,end) & ~active);
   B       = [descent(:,Bgood), I(:,Bbad), -I(:,Bbad)];
   DN      = I(:,Nbad);
   dhat    = -sum([-descent(:,Ngood), DN],2);
else
   [B,N]   = getTangentCone(iterate.x,O,tol,scale,scheme,success);
   Bgood   = (B'*infGrad(:,end) == 0 & B'*g < 0);
   Ngood   = (N'*infGrad(:,end) == 0);
   Bbad    = setdiff(1:size(B,2),Bgood);
   Nbad    = setdiff(1:size(N,2),Ngood);
   B       = [B(:,Bgood), B(:,Bbad), -B(:,Bbad)];
   DN      = N(:,Nbad);
   NN      = [N, -N];
   dhat    = -sum([NN(:,NN'*g > 0), DN],2);
end
if ~any(dhat) || isequal(d(:,end),dhat), dhat = []; end
if ~strcmp(strategy,'Gradient_3n2n'), descent = []; end
D = [scaleDirections(d,scale), dhat, scaleDirections(descent,scale),B,DN];
D(:,~any(D)) = [];
return

%*******************************************************************************
% getTangentCone:  Compute and return tangent cone generators.
% ------------------------------------------------------------------------------
% Acknowledgement:  The mathematics associated with the more basic parts of this 
%    algorithm are due to R. M. Lewis & V. Torczon (College of William & Mary)
% ------------------------------------------------------------------------------
% Called by: standDirections, gradDirections
% Calls:     activeConstraints, removeRedundancy,
%            isFullRank (in-line), isDegenerate (in-line)
% VARIABLES:
%  B        = V*inv(V'*V)
%  N        = vectors than span the null space of V'
%  x        = current iterate continuous variables
%  Omega    = structure defining feasible region Omega = {x: l<=A*x<=u}
%    .A     =   matrix of linear constraint coefficients
%    .l     =   vector of constraint lower bounds
%    .u     =   vector of constraint upper bounds
%  tol      = tolerance within which constraints are considered active
%  scale    = scale factors for scaling the appropriate Poll directions
%  scheme   = scheme for handling degenerate nonredundant linear constraints
%  success  = flag indicating previous iteration was successful
%  n        = number of continuous variables
%  I        = identity matrix
%  Ax       = A*x
%  tolFinal = final tolerance to determine e-active constraints
%  V        = matrix whose columns are normal to active linear constraints
%  lactive  = indices of linear constaints active at their lower bound
%  uactive  = indices of linear constaints active at their upper bound
%  active   = indices of reformulated linear constraints
%  A        = alternative storage of linear constraint coefficients
%  b        = alternative storage of linear constraint bounds
%  r        = rank of A after linearly dependent rows are deleted
%  m        = number of rows in A after linearly dependent rows are deleted
%  ind      = indices of constraints used to form tangent cone generators
%  AA       = A(ind,:)
%  [Q,R]    = QR factorization of V
%*******************************************************************************
function [B,N] = getTangentCone(x,Omega,tol,scale,scheme,success)

isDegenerate = @(A) rank(A) < size(A,1);
isFullRank   = @(A) rank(A) == min(size(A));

% Initialize variables
n  = length(x);
I  = eye(n);
Ax = Omega.A*x;
tolFinal = tol^2/2;

% Check for degeneracy, reducing the tolerance as a remedy
V = 0;
while isDegenerate(V')
   lactive = activeConstraints( Ax,  Omega.l,  Omega.u, tol);
   uactive = activeConstraints(-Ax, -Omega.u, -Omega.l, tol);
   V = [ Omega.A(find(uactive),:)', -Omega.A(find(lactive),:)' ];  %#ok
   if tol <= tolFinal, break, end
   tol = tol/2;
end

% If degenerate, remove any redundant constraints
if isDegenerate(V')
   A = [Omega.A; -Omega.A]; 
   b = [Omega.u; -Omega.l];
   active = activeConstraints(A*x, [Omega.l; -Omega.u], b, tol);
   [A,b]  = removeRedundancy(A, b, x, active);

   % If still degenerate, compute V for a different nearby cone
   m = length(b);
   r = rank(A);
   if r < m
      normsA = vecnorm(A,2,2);
      dist = (b - A*x)./normsA;
      V = 0;
      while ~isFullRank(V)
         if strcmp(scheme,'full')
            done = 0;
            B = []; N = [];
            while ~done
               ind = getAlternateCone('sequential',m,r,dist,0);
               if isequal(ind,1:r) && ~isempty(N)
                  done = 1;
               else
                  V = A(ind,:)';
                  if isFullRank(V)
                     [Q,R] = qr(V,0);
                     tempB = Q/R';
                     tempN = I - tempB*V';
                     B = [B, tempB];  %#ok
                     N = [N, tempN];  %#ok
                  end
               end
            end
            N(abs(N) < 10*eps) = 0;   %#ok
            B(abs(B) < 10*eps) = 0;   %#ok
            return
         else
            ind = getAlternateCone(scheme,m,r,dist,success);
            V = A(ind,:)';
         end
      end
   end
end

% If still degenerate, there must be an error somewhere
if ~isFullRank(V)
   error('mads:constraints:degenerate', ...
         'Unresolved degenerate linear constraints (getTangentCone).');
end

% Compute B and N
if isempty(V)
   B = zeros(n,1);
   B(:,1:end) = [];
   N = scaleDirections(I,scale);
else
   [Q,R] = qr(V,0);
   B = Q/R';
   N = I - B*V';
end

% Zero out elements that are very small already
N(abs(N) < 10*eps) = 0;
B(abs(B) < 10*eps) = 0;
return

%*******************************************************************************
% getAlternateCone:  Gets generators of alternative cone if degeneracy occurs.
% ------------------------------------------------------------------------------
% Called by: getTangentCone
% VARIABLES:
%  ind     = indices of constraints to be used in computing tangent cone
%  scheme  = method used to choose constraints
%  m       = number of e-active, nonredundant linear constraints
%  r       = rank of A
%  dist    = vector of distances from current point to constraint boundaries
%  success = flag indicating previous iteration was successful
%  extra   = an extra index added randomly
%  y       = temporary storage
%*******************************************************************************
function ind = getAlternateCone(scheme,m,r,dist,success)

switch scheme

   % Select indices based on distance to constraint boundary
   case 'closest'
      [~,ind] = sort(dist);
      extra   = randperm(m-r);
      ind     = [ind(1:r-1)', ind(r-1+extra(1))'];

   % Random Selection of constraint indices
   case 'random'
      ind = randperm(m);
      ind = ind(1:r);

   % Sequential Selection of constraint indices
   case 'sequential'
      ind = getappdata(0,'ALTCONE');
      if length(ind) ~= r || max(ind) > m || isequal(ind,(m-r+1:m))
         ind = 1:r;
      elseif success == -1
         for i = 1:r
            if ind(r+1-i) < m+1-i
               k = ind(r+1-i);
               for j = r+1-i:r
                  ind(j) = k+j+i-r;
               end
               break
            end
         end
      end
      setappdata(0,'ALTCONE',ind);

   % Error if the scheme is not from the list
   otherwise
      error('mads:constraints:degenerate', ...
            'Invalid scheme for choosing constraints (getAlternativeCone).');
end
return

%*******************************************************************************
% scaleDirections:  Scale each column by a vector of scale factors.
% ------------------------------------------------------------------------------
% Called by: standDirections, madsDirections, gradDirections
% VARIABLES:
%  Y     = output matrix
%  X     = input matrix
%  scale = vector of scale factors
%*******************************************************************************
function Y = scaleDirections(X,scale)
if isempty(X) || all(scale == 1.0)
   Y = X;
else
   Y = diag(scale)*X;
end
return

%*******************************************************************************
% activeConstraints:  Check which constraint bounds are numerically active.
% ------------------------------------------------------------------------------
% Called by: standDirections, gradDirections, getTangentCone
% VARIABLES:
%  active = logicals indicating which constraints is numerically active
%  x      = point in R^n to be tested 
%  [a,b]  = upper and lower bounds on x
%  tol    = error tolerance for e-active constraints
%  infa   = indices of x having no lower bound
%  infb   = indices of x having no upper bound, but having a lower bound
%  finite = indices of x having finite lower and uppr bounds
%*******************************************************************************
function active = activeConstraints(x,a,b,tol)

% Set status of each continuous variable
% infa   = find(isinf(a));
infb   = find(isinf(b) & ~isinf(a));
finite = isfinite(a.*b);
one    = ones(length(infb),min(1,length(infb)));

% Determine status of each variable
active         = zeros(length(x),1);
scale          = max([abs(x(infb)),abs(a(infb)),one],[],2);
active(infb)   = (x(infb)  -a(infb)   <= tol*scale);
scale          = max(b(finite)-a(finite),1);
active(finite) = (x(finite)-a(finite) <= tol*scale);
return

%*******************************************************************************
% removeRedundancy:  Removes rows of Ax <= b that are redundant.
% ------------------------------------------------------------------------------
% This algorithm is courtesy of Olga A. Brezhneva (Miami University, Ohio)
% ------------------------------------------------------------------------------
% Called by: getTangentCone
% Calls:     runLP
% VARIABLES:
%  AA     = matrix A with only nonredundant e-active constraints left
%  bb     = vector b with only nonredundant e-active constraints left
%  A      = coefficient matrix
%  b      = right-hand side
%  x      = current iterate continuous variables
%  active = indices of the active constraints (rows) in A
%  ind    = vector of indices
%  bb     = b vector with only nonredundant e-active constraints left
%  m      = length of bb
%  res    = residual of current constraint
%  px     = projection of x onto the current constraint boundary
%  pres   = vector of residuals with respect to px
%  NR     = logical indicating non-redundancy status
%  xval   = optimal x with respect to the LP solver
%  fval   = LP objective function value at xval
%  i      = constraint counter
%*******************************************************************************
function [AA,bb] = removeRedundancy(A,b,x,active)

% Form e-active constraints
ind = find(active);
AA  = A(ind,:);
bb  = b(ind);

% Begin loop through constraints
i = 1;
m = length(bb);
while i <= m
   
   % Strategy 1: verification of nonredundancy
   res  = bb(i) - AA(i,:)*x;
   px   = x + AA(i,:)'*res;
   pres = bb - AA*px;
   pres(i) = 1;
   NR = all(pres > 0);

   % Strategy 2: solve LP problem
   if ~NR
      ind    = 1:length(bb);
      ind(i) = [];
      [~,fval] = runLP(-AA(i,:)',AA(ind,:),bb(ind));
      if isempty(fval)
         NR = 1;
      else
         NR = (-fval > bb(i));
      end
   end

   % Delete constraint if redundant
   if NR
      i = i + 1;
   else
      AA(i,:) = [];
      bb(i)   = [];
      m = length(bb);
   end
end
return

%*******************************************************************************
% dVector:  Returns the appropriate 3^n vector approximating v.
% ------------------------------------------------------------------------------
% Called by: gradDirections
% VARIABLES:
%  w = L-1, closest L-2, or L-Inf {-1,0,1}^n approximation to v
%  v = input vector
%  p = norm of the gradient
%*******************************************************************************
function w = dVector(v,p)
switch p
case {1}
   w = (abs(v) == max(abs(v))).*sign(v);
case {2}
   w = v/norm(v,inf);
   w(~isfinite(w) | abs(w) < tan(pi/8)) = 0;
   w = sign(w);
case {inf}
   w = sign(v);
end
return

%*******************************************************************************
% getPollOrder:  Returns the order in which to poll from the poll set.
% ------------------------------------------------------------------------------
% Called by: poll
% VARIABLES:
%  order     = resulting Poll order
%  orderCode = string indicating user choice
%  D         = set of Poll directions
%  P         = the Poll set
%  RunData   = structure of Run statistics and information
%    .goodD  =   previous direction of the last successful iteration
%    .porder =   poll order from the previous successful iteration
%  n         = dimension of the optimization problem
%  nD        = number of Poll directions
%  angles    = a measure of angles between .goodD or .sd and the columns of D
%  sf,sc,sh  = surrogate f-, c-, and h-values
%*******************************************************************************
function order = getPollOrder(orderCode,D,P,RunData,surrogate)

% get the dimensions of D
[n,nD] = size(D);
order  = 1:nD;

% Construct the order according to user choice
switch orderCode
   case 'Alternating'
   switch nD
      case {n+1}
         order = [nD, 1:nD-1];
      case {2*n}
         order = reshape([1:n; n+1:2*n],1,nD);
      otherwise
         order = 1:nD;
   end
   case 'Random'
      order = randperm(nD);
   case 'Dynamic'
      if isfield(RunData,'porder')
         angles  = (RunData.goodD'*D) ./ sqrt(sum(D.^2));
         [~,ind] = min(-angles);
         k       = find(RunData.porder == ind);
         order   = [ind, RunData.porder(1:k-1), RunData.porder(k+1:end)];
      end
   case 'DynamicRanked'
      if isfield(RunData,'porder')
         angles    = (RunData.goodD'*D) ./ sqrt(sum(D.^2));
         [~,order] = sort(-angles);
      end
   case 'SimplexGradient'
      if isfield(RunData,'sd') && ~isempty(RunData.sd)
         angles    = (RunData.sd'*D) ./ sqrt(sum(D.^2));
         [~,order] = sort(-angles);
      end
   case 'Surrogate'
      sf = feval(surrogate.evaluator,[P.x]',surrogate.f);
      if isfield(surrogate,'c')
         nC = length(surrogate.c);
         sc = zeros(n,nC);
         for k = 1:nC
            sc(:,k) = feval(surrogate.evaluator,[P.x]',surrogate.c(k));
         end
         sh = sum((sc > 0).*sc.*sc,2);
      end
      if RunData.pcenter.h > 0
         [~,order] = sort(sh);
      else
         [~,order] = sort(sf);
      end
   case 'Custom'
      order = RunData.porder;
      if length(order) ~= nD
         error('mads:poll:dim', ...
              ['Poll directions and orders have incompatible dimensions ', ...
               '(getPollOrder).']);
      end
end
return

%*******************************************************************************
% FUNCTIONS CALLED BY MADS UPDATE ROUTINE
%*******************************************************************************
%*******************************************************************************
% getPollCenter: Gets the POLL center based on CenterCode.
% ------------------------------------------------------------------------------
% Called by: processInput, processOutput, search, mvpPoll, update
% VARIABLES:
%  pcenter      = chosen POLL center
%  iCache       = Cache index of poll center
%  centerCode   = code indicating which points to POLL around
%                 (0 = current solution; otherwise = Filter(centerCode))
%  Filter       = sorted indices of Cache iterates in the filter
%  Cache        = collection of all previous iterates
%  closest      = min(centerCode, filter size)
%  noFilter     = flag indicating no infeasible iterates have been found
%  noFeasible   = flag indicating no feasible iterates have been found
%  filterPoints = indices of valid filter points
%*******************************************************************************
function [pcenter,iCache] = getPollCenter(centerCode,Filter,Cache,rho)

% Error check existence of a possible poll center.
noFilter   = isempty(Filter.F);
noFeasible = isempty(Filter.feasible);
if noFeasible && noFilter
   error('mads:pollCenter:badFilter','No suitable poll center (PollCenter).');
end

% Progressive Barrier Poll Center
if isempty(centerCode)   
   [pcenter,iCache] = getProgBarrierPollCenter(Filter,Cache,rho);
   return
end

% Delete repeated poll center codes, and retrieve infeasible poll center codes
centerCode  = unique(centerCode);
filterCode  = nonzeros(centerCode);
nFilterCode = length(filterCode);

% Process any request for best feasible point
iCache = [];
if nFilterCode < length(centerCode)
   if noFeasible
      if ~(any(filterCode) == 1)
         iCache = Filter.LIP;
      end
   else
      iCache = Filter.BFP;
   end
end

% Process requests for filter points
if ~isempty(filterCode)
   if noFilter
      if isempty(iCache)
         iCache = [iCache, Filter.BFP];
      end
   else
      filterPoints = find([Cache.iterate(Filter.F).h] < Filter.hmax);
      if isempty(filterPoints)
         closest = 1;
      else
         closest = min([filterCode(:)'; ...
                        filterPoints(end)*ones(1,nFilterCode)],[],1);
      end
      iCache = [iCache, Filter.F(unique(closest))];
   end
end

% Want particular filter point.  If none, pick current best feasible point.
% for k = 1:length(filterCode)
%    if noFilter
%       iCache = Filter.feasible(1);
%    else
%       filterPoints = find([Cache.iterate(Filter.F).h] < Filter.hmax);
%       if isempty(filterPoints)
%          closest = 1;
%       else
%          closest = min(filterCode(k),filterPoints(end));
%       end
%       iCache = Filter.F(closest);
%    end
% end
pcenter = Cache.iterate(iCache);
return

%*******************************************************************************
% getProgBarrierPollCenter: Gets the POLL center based the progressive barrier.
% ------------------------------------------------------------------------------
% Called by: getPollCenter
% VARIABLES:
%  pcenter      = chosen POLL center
%  iCache       = Cache index of poll center
%  centerCode   = code indicating which points to POLL around
%                 (0 = current solution; otherwise = Filter(centerCode))
%  Filter       = sorted indices of Cache iterates in the filter
%  Cache        = collection of all previous iterates
%  closest      = min(centerCode, filter size)
%  noFilter     = flag indicating no infeasible iterates have been found
%  noFeasible   = flag indicating no feasible iterates have been found
%  filterPoints = indices of valid filter points
%*******************************************************************************
function [pcenter,iCache] = getProgBarrierPollCenter(Filter,Cache,rho)

if isempty(Filter.IIP)
   iCache = Filter.BFP;
elseif isempty(Filter.BFP)
   iCache = Filter.IIP;
elseif Cache.iterate(Filter.BFP).f - rho > Cache.iterate(Filter.IIP).f
   iCache  = [Filter.IIP, Filter.BFP];
else
   iCache  = [Filter.BFP, Filter.IIP];
end
pcenter = Cache.iterate(iCache);

return


%*******************************************************************************
% updateOmega:  Update the set of linear constraints .l <= .A*x <= .u.
% ------------------------------------------------------------------------------
% Parts of this routine are courtesy of Dennis & Schnabel textbook
% ------------------------------------------------------------------------------
% Called by: processInput, search, mvpPoll, update, evalPointSet
% Calls:     getpID, < Omega File >
% VARIABLES:
%  Omega               = structure describing (linear) feasible region
%    .A                =   coefficient matrix
%    .l                =   vector of lower bounds
%    .u                =   vector of upper bounds
%    .plist            =   list of allowable categorical variable values
%    .s                =   scale factors for variable scaling
%    .r                =   shifts for variable scaling
%  newScale            = scale factors for constructing mesh
%  Problem             = structure of problem data
%    .File.O           =   name of file defining feasible region Omega
%    .isMVP            =   flag indicating problem is an MVP problem
%  Options             = structure of MADS options
%    .removeRedundancy =   turns on/off removal of redundant constraints
%    .scale            =   base for logarithmic scaling
%  iterate             = current iterate
%    .n                =   dimension of the continuous variable space
%    .x                =   vector of continuous variable values
%    .p                =   cell array of categorical variables
%  n                   = number of continuous variables of iterate
%  m                   = number of linear constraints
%  [a,b]               = temporary storage of variable bounds
%  scaleA              = vector of scales for scaling linear constraints
%  AA,bb               = reformulation of linear constraints
%  Redundant           = flags indicating redundant linear constraints
%  ind                 = indices of nonredundant linear constraints
%  indcopy             = temporary copy of ind
%  xval,fval           = x-value and f-value of LP solution
%**************************************************************************
function [Omega,newScale] = updateOmega(Problem,Options,iterate)

% Get parameters that define the linear constraints
n = length(iterate.x);
if exist(Problem.File.O,'file')
   switch nargin(Problem.File.O) + 2*Problem.isMVP
   case {1}
      [A,l,u] = feval(Problem.File.O,n);
   case {2}
      error('mads:constraints:badFileFormat', ...
            'Non-MVP Omega file cannot have a second argument (updateOmega).');
   case {3}
      [A,l,u,plist] = feval(Problem.File.O,n);
   case {4}
      [A,l,u,plist] = feval(Problem.File.O,n,iterate.p);
   end
else
   if Problem.isMVP
      error('mads:constraints:file', ...
           ['MVP problem missing its Omega file, ',Problem.File.O, ...
            ' (updateOmega).']);
   else
      [A,l,u] = deal(eye(n),-Inf*ones(n,1),Inf*ones(n,1));
      if isfield(Problem,'Param') && isfield(Problem.Param,'A')
         A = Problem.Param.A;
      end
      if isfield(Problem,'Param') && isfield(Problem.Param,'l')
         l = Problem.Param.l;
      end
      if isfield(Problem,'Param') && isfield(Problem.Param,'u')
         u = Problem.Param.u;
      end
   end
end
[m,n] = size(A);

if ~isequal(size(l),size(u),[m,1]) || any(l > u) 
   error('mads:constraints:dim','Invalid constraint dimensions (updateOmega).');
end

% Get shifts for variable scaling
if Options.scaleVariables
   [s,r]   = getScaleShift(l(1:n),u(1:n));
   [A,l,u] = scaleOmega(A,l,u,s,r);
   Options.scale = 0;
else
   s = ones(n,1);
   r = zeros(n,1);
end

% Scale the linear constraints by their inf-norms and form Omega structure
a = l(1:n);
b = u(1:n);
scaleA = max(abs(A),[],2);
A = diag(1./scaleA)*A;
l = l./scaleA;
u = u./scaleA;

% Detect redundant linear constraints (due in part to Olga Brezhneva)
if Options.removeRedundancy
   AA = [-A; A];
   bb = [-l; u];
   Redundant = ~isfinite(bb)';
   ind = find(isfinite(bb));
   for i = 1:length(ind)
      indcopy = ind;
      indcopy(i) = [];
      [~,fval] = runLP(-AA(ind(i),:)',AA(indcopy,:),bb(indcopy));
      if isempty(fval)
         Redundant(ind(i)) = 0;
      else
         Redundant(ind(i)) = (-fval <= bb(ind(i)));
      end
   end

   % Delete redundant linear constraints
   ind = all(reshape(Redundant,m,2)');
   A(ind,:)   = [];
   l(ind)     = [];
   u(ind)     = [];
end

% Form Omega structure and get new scale factors
Omega = struct('A',full(A),'l',l,'u',u,'s',s,'r',r);
if Problem.isMVP
   Omega.plist = plist;
   newScale = getScaleFactors(Options.scale,iterate,a,b,Problem.nameCache, ...
                              plist);
else
   newScale = getScaleFactors(Options.scale,iterate,a,b,Problem.nameCache);
end
return

%*******************************************************************************
% getScaleShift:  Compute variable scaling and shift parameters. 
% ------------------------------------------------------------------------------
% Called by: updateOmega
% VARIABLES:
%  s,r = scale factors and shifts for variable scaling option
%  l,u = variable lower and upper bounds
%  n   = dimension of feasible region
%*******************************************************************************
function [s,r] = getScaleShift(l,u)
n = length(u);
r = zeros(n,1);
s = ones(n,1);
ind = l > -inf;
r(ind) = l(ind);
ind = isfinite(u.*l);
r(l > -inf) = l;
s(ind) = u(ind) - l(ind);
return

%*******************************************************************************
% scaleOmega:  Compute variable scaling and shift parameters. 
% ------------------------------------------------------------------------------
% Called by: updateOmega
% VARIABLES:
%  s,r = scale factors and shifts for variable scaling option
%  l,u = variable lower and upper bounds
%  n   = dimension of feasible region
%*******************************************************************************
function [newA,newl,newu] = scaleOmega(A,l,u,s,r)
%n = size(A,2);
newA = A*diag(s);
newl = l - A*r;
newu = u - A*r;
return

%*******************************************************************************
% getScaleFactors:  Get scale factors for scaling directions.
% ------------------------------------------------------------------------------
% Called by: UpdateOmega
% VARIABLES:
%  newScale  = vector of scale factors for scaling of directions
%  scale     = base of logarithmic scaling (0 = no scale)
%  iterate   = current poll center
%  [a,b]     = vectors of lower and upper bounds on continuous variables
%  isMVP     = flag indicating problem is an MVP problem
%  nameCache = name of the base workspace Cache variable
%  typx      = "typical" values used for scaling if bounds are insufficient
%  pID       = ID for the categorical variable values
%  ind       = index used to match Cache points with the same pID
%  zeroInf   = flags for zero or infinite values in [a,b]
%  bad       = indices of variables with a and b both zero or infinite
%  good      = indices of variables with a and b both finite and nonzero
%  lhalf     = indices of variables with finite nonzero a and zero or infinite b
%  rhalf     = indices of variables with finite nonzero b and zero or infinite a
%  expScale  = exponents of the approximate ranges of the variables
%  baseScale = factor for coverting from base b to base 10 logarithm
%*******************************************************************************
function newScale = getScaleFactors(scale,iterate,a,b,nameCache,varargin)

% Compute scale factors, if selected
n = length(iterate.x);
if scale

   % Compute "typical" value of x for scaling
   isMVP = (nargin > 5);
   Cache = getappdata(0,nameCache);
   if isempty(Cache.iterate)
      typx = iterate.x;
   else
      if isMVP
         pID = getpID(iterate.p,varargin{:});
         ind = find(Cache.pID(1:Cache.size) == pID);
         if isempty(ind)
            typx = iterate.x;
         else
            typx = Cache.iterate(ind(1)).x;
         end
      else
         typx = Cache.iterate(1).x;
      end
   end

   % Flag zero or infinite bounds, and categorize each variable's bounds
   zeroInf = [ (~isfinite(a)|abs(a)<=eps), (~isfinite(b)|abs(b)<=0)];
   good    = (~any(zeroInf,2));
   lhalf   = (~zeroInf(:,1)  &  zeroInf(:,2));
   rhalf   = ( zeroInf(:,1)  & ~zeroInf(:,2));
   bad     = (all(zeroInf,2) &  typx);
   verybad = (all(zeroInf,2) & ~typx);

   % Scale variables with finite nonzero bounds IAW D&S, ch 7 and pp. 278-9
   baseScale       = 1/log10(scale);
   expScale        = ones(n,1);
   expScale(good)  = baseScale * (log10(abs(a(good)))+log10(abs(b(good))))/2;
   expScale(lhalf) = baseScale *  log10(abs(a(lhalf)));
   expScale(rhalf) = baseScale *  log10(abs(b(rhalf)));
   expScale(bad)   = baseScale *  log10(abs(typx(bad)));
   expScale(verybad) = 0;
   newScale = scale.^round(expScale);
else
   newScale = ones(n,1);
end
return

%*******************************************************************************
% getSimplexGradient:  Get simplex gradient at current iterate.
% ------------------------------------------------------------------------------
% Called by: update
% VARIABLES:
%  sd         = simplex gradient at the current iterate
%  lambda     = lambda used in lambda-poisedness check
%  Cache      = structure of previously computed iterates
%    .iterate = vector of previously computed iterates
%      .x     =   continuous variable values
%      .f     =   objective function value
%  iCenter    = index into Cache corresponding to current poll center
%  delta      = current poll size
%  success    = 0-1 indicator of success of the previous iteration
%  ind        = vector indices that get redefined in every block of this code
%  X          = matrix storage of continuous variable values of all iterates
%  S          = system matrix used in computing simplex derivatives
%*******************************************************************************
function sg = getSimplexGradient(lambda,Cache,iCenter,delta,success)

sg = [];

% Remove poll center from set of cache points
ind           = [1:iCenter-1, iCenter+1:Cache.size];
pcenter       = Cache.iterate(iCenter);
pID           = Cache.pID(iCenter);
Cache.pID     = Cache.pID(ind);
Cache.iterate = Cache.iterate(ind);

% Filter cache if MVP problem, and return if there are not enough data points
ind     = (Cache.pID == pID);
Cache   = Cache.iterate(ind);
n       = length(pcenter.x);
if length(Cache) < n + 1, return, end

% Translate points to unit ball centered at current poll center
X = [Cache(2:end).x];
F = [Cache(2:end).f];
S = (X - diag(pcenter.x)*ones(size(X)))/(2*delta);

% Filter out all points not in the trust region
ind     = find(sum(S.^2) <= 1);
nPoints = length(ind);
if nPoints < n, return, end
S = S(:,ind);
F = F(ind);

% Check for block lambda-poisedness
ind    = nPoints-n+1:nPoints;
Y      = S(:,ind)';
poised = cond(Y)/norm(Y) <= lambda;

% Check for lambda-poisedness one-by-one
if ~poised
   ind = [];
   for k = nPoints:-1:1
      Y  = S(:,[k,ind])';
      poised = cond(Y)/norm(Y) <= lambda;
      if poised
         ind = [k, ind];  %#ok
      end
   end
   if length(ind) < n, return, end
end

% Construct simplex gradient from a set of lambda-poised points
A = S(:,ind)';
b = F(ind)' - pcenter.f;
if size(A,1) < size(A,2)
   sg = lsqminnorm(A,b);
else
   sg = A\b;
end
return


%*******************************************************************************
% FUNCTIONS FOR EVALUATING FUNCTIONS AND UPDATING CACHE DATA
%*******************************************************************************
%*******************************************************************************
% evalPointSet:  Evaluates the objective and constraint functions.
% ------------------------------------------------------------------------------
% Called by: processInput, search, poll, mvpPoll
% Calls:     updateOmega, inOmega, getpID, isCacheHit, evalFunction,
%            updateCache, updateFilter
% VARIABLES:
%  P              = set of iterates to be evaluated
%  success        = logical indicating a successful iteration
%  Filter         = updated filter
%  ptype          = type of function evaluation (e.g., S=search,P=poll, etc.)
%  Problem        = structure containing problem parameters
%    .File.O      =   name of Omega file
%    .fType       =   type of functions file (C=C, F=FORTRAN, M=MATLAB)
%    .isMVP       =   logical indicating if problem is n MVP
%    .Omega       =   structure defining feasible region {x: l<=Ax<=u}
%      .A         =     coefficient matrix
%      .l         =     vector of lower bounds
%      .u         =     vector of upper bounds
%      .plist     =     list of valid categorical variables
%    .nameCache   =   name of the base workspace Cache variable
%  full           = flag indicating if all points in P will be evaluated
%  Options        = structure of user options
%    .computeGrad =   flag for computing gradients
%  Filter         = index into the Cache of (nondominated) points in the filter
%  nP             = number of iterates to be evaluated
%  pID            = categorical variable value ID
%  xnorm          = 1-norm of iterate being evaluated
%  hit            = flag indicating a Cache hit
%  unfiltered     = logical indicating an iterate is unfiltered
%*******************************************************************************
function [P,success,Filter] = evalPointSet(ptype,Problem,P,full,Options,...
                                              Filter,add2sur)
success = -1;
if isempty(P), return; end
Cache  = getappdata(0,Problem.nameCache);
nCache = Cache.size;
grad   = Options.computeGrad;

% Initialize output fields
nP = length(P);
for k = 1:nP
   P(k).n    = length(P(k).x);
   P(k).type = ptype;
end
[P.sf,P.sc,P.sh,P.f,P.c,P.h,P.param] = deal([]);

% Initialize output fields used with derivatives
if grad
   [P.sgradf,P.sgradc,P.sgradh,P.gradf,P.gradc,P.gradh] = deal([]);
end

% Evaluate each iterate, as appropriate
for k = 1:nP

   % update MVP Omega parameters, as is necessary
   if Problem.adaptOmega
      [Problem.Omega] = updateOmega(Problem,Options,P(k));
   end

   % Store variable values in strings for debugging
   if Options.debug >= 3
      disp('   ---------------------------');
      disp(['   evalPointSet:  type = ',ptype]);
   end
   if ~isempty(P(k).p)
      p = '';
      if Options.debug >= 4, disp('   evalPointSet:  p = '); end
      for i = 1:length(P(k).p)
         if ischar(P(k).p{i})
            p{i} = P(k).p{i};
         else
            p{i} = int2str(P(k).p{i});
         end
          if Options.debug >= 4, disp(['     ',[p{i}]]); end
      end
   end
   if Options.debug >= 4
      disp('   evalPointSet:  x = ');
      disp([repmat('     ',P(k).n,1),num2str(P(k).x,Options.digits)]);
   end

   % Test if iterate lies in Omega
   P(k).f = NaN;
   if inOmega(P(k),Problem.Omega)

      % Give trial point unique IDs to expedite Cache searching
      pID = 0;
      if Problem.isMVP
         pID = getpID(P(k).p,Problem.Omega.plist);
      end
      xnorm = norm(P(k).x,inf);

      % Test to see if the iterate is in the Cache
      Cache       = getappdata(0,Problem.nameCache);
      hit         = isCacheHit(Cache,P(k),xnorm,pID,50000);
      Cache.nHits = Cache.nHits + hit;
      setappdata(0,Problem.nameCache,Cache);

      % Evaluate function, store in Cache, and update Cache and filter
      if hit
         if Options.debug >= 3, disp('   evalPointSet: *** Cache hit ***'); end
      else
         
         % Grab constant categorical variable values if that option was chosen
         if Options.fixCategorical
            P(k).p = Problem.Param.p;
         end

         % Closed constraint evaluation (X-feasible check)
         if isempty(Problem.File.X)
            feasible = 1;
         else
            [P(k),feasible] = evalClosedConstraints(Problem.File.X,P(k),grad);
            if Options.debug >= 3
               disp(['   evalPointSet:  X-feasible = ',int2str(feasible)]);
            end
         end

         % Evaluate objective/constraints if initial points or X-feasible
         if feasible || strcmp(ptype,'0')
            P(k) = evalFunction(Problem.File.F,Problem.fType,P(k),grad);
            if Options.debug >= 3
               disp(['   evalPointSet:  f = ',num2str(P(k).f,Options.digits)]);
            end
         end

         % Update the cache
         if Options.fixCategorical
            P(k).p = {};
         end
         if feasible || strcmp(ptype,'0') || Options.countNotInX
            Cache = updateCache(Cache,Problem,P(k),xnorm,pID,1,0);
            setappdata(0,Problem.nameCache,Cache);
         end

         % Update filter
         if isfinite(P(k).h)
            [unfiltered,Filter] = updateFilter(P(k),Filter,Cache);
         else
            unfiltered = 0;
         end
         if unfiltered
            success = 1;
         end
         if add2sur
            Cache = getappdata(0,Problem.nameCache);
            Cache.isSurrPoint(Cache.size) = 1;
            setappdata(0,Problem.nameCache,Cache);
         end
         if success >= 0 && ~full, break; end
      end
   else
      if Options.debug >= 3, disp('   evalPointSet: *** not in Omega ***'); end
   end
   if Options.debug >= 3, disp('   ---------------------------'); end
end
Cache = getappdata(0,Problem.nameCache);
if Options.debug >= 2 
   disp(['   Points identified/evaluated: ', ... 
                           int2str(nP),' / ',int2str(Cache.size-nCache)]);
end
return

%*******************************************************************************
% getpID:  Get unique ID based on values of categorical variables.
% ------------------------------------------------------------------------------
% Called by: evalPointSet
% VARIABLES:
%  pID   = unique ID
%  p     = cell array of categorical variables
%  plist = cell array of allowable values for each categorical variable
%  m     = vector of plist lengths for each categorical variable
%  ind   = index of plist that matches the current categorical variable
%*******************************************************************************
function pID = getpID(p,plist)

np  = length(p);
m   = zeros(np,1);
pID = zeros(np,1);
for k = 1:np
   m(k) = length(plist{k});
   if ischar(p{k})
      ind = find(strcmp(p{k},plist{k}));
   else
      ind = find([plist{k}{:}] == p{k});
   end
   pID(k) = ind(1);
end
pID = pID'*[1; cumprod(m(1:end-1)+1)];
return

%*******************************************************************************
% inOmega:  Test to see if point x satisfies l <= A*x <=u.
% ------------------------------------------------------------------------------
% Called by: makeFeasible, evalPointSet
% VARIABLES:
%  pass     = logical indicating is x satisfies l <= A*x <=u
%  iterate  = current iterate
%    .x     =   vector of continuous variables
%    .p     =   cell array of categorical variables
%  Omega    = feasible space with respect to bound and linear constraints
%    .A     = coefficient matrix for linear and bound constraints
%    .l     = lower bound vector for linear and bound constraints
%    .u     = upper bound vector for linear and bound constraints
%    .plist = cell array of lists of allowed categorical values in p
%    .s     = variable scale factors
%    .r     = variable shift parameters
%  Ax       = .A*x
%*******************************************************************************
function pass = inOmega(iterate,Omega)

% Scale variables to adjusted linear constraints
x = iterate.x./Omega.s - Omega.r;

% Check for errors in categorical variables or Omega construction
for k = 1:length(iterate.p)
   if ischar(iterate.p{k})
      pass = any(strcmp(iterate.p{k},Omega.plist{k}));
   else
      pass = any([Omega.plist{k}{:}] == iterate.p{k});
   end
   if ~pass, return; end
end

% Check if iterate is in Omega
Ax   = Omega.A*x;
pass = all(Omega.l <= Ax) && all(Ax <= Omega.u);
return

%*******************************************************************************
% isCacheHit:  Test to see if point has already been evaluated previously.
% ------------------------------------------------------------------------------
% Called by: evalPointSet
% VARIABLES:
%  hit        = logical indicating that the iterate matches a Cache point
%  Cache      = array of previously computed iterates
%    .iterate =   vector of iterates
%    .tol     =   tolerance for declaring a Cache hit w.r.t iterate.x
%    .xnorm   =   vector of 1-norms of iterates in the Cache
%    .pID     =   vector of IDs of categorical variable values
%  iterate    = current point
%    .x       =   coordinates of continuous  variables
%    .p       =   coordinates of categorical variables
%  xnorm      = 1-norm of iterate.x
%  pID        = ID of categorical variable values
%  looksize   = how far back into the Cache will be searched
%  ind        = indices of Cache points that may match the current point
%  xvalues    = temporary storage of Cache continuous variable values
%*******************************************************************************
function hit = isCacheHit(Cache,iterate,xnorm,pID,looksize)

% Initialize Cache search
hit = 0;
ind = max(1,Cache.size-looksize):Cache.size;

% Screen out iterates without sufficiently close x-norms (see Apostol, p49)
ind = ind((abs(Cache.xnorm(ind) - xnorm) < Cache.tol));
if isempty(ind), return, end

% Screen out iterates with different categorical variable values
if pID ~= 0
   ind = ind(Cache.pID(ind) == pID);
   if isempty(ind), return, end
end

% Finish search
xvalues = [Cache.iterate(ind).x];
ind = ind(max(abs(xvalues - repmat(iterate.x,1,length(ind))),[],1) < Cache.tol);
hit = ~isempty(ind);
return

%*******************************************************************************
% evalClosedConstraints:  Evaluate closed constraints file
% ------------------------------------------------------------------------------
% Called by: evalPointSet
% Calls: < User X >
% VARIABLES:
%  iterate  = current iterate
%    .x     =   coordinates of iterate
%    .p     =   values of categorical variables
%    .f     =   objective function value of iterate
%    .c     =   constraint function values of iterate
%    .h     =   constraint violation function value of iterate
%    .gradf =   f'(x) at current iterate x
%    .gradc =   c'(x) at current iterate x
%    .gradh =   h'(x) at current iterate x
%  feasible = logical indicating if closed constraints are satisfied
%  xFile    = name of closed constraints file, if it exists
%  grad     = flag for computing gradients, if available
%*******************************************************************************
function [iterate,feasible] = evalClosedConstraints(xFile,iterate,grad)

switch nargin(xFile)
   case {1}
      feasible = feval(xFile,iterate.x);
   case {2}
      feasible = feval(xFile,iterate.x,iterate.p);
   otherwise
      error('mads:function:NumberOfArguments',...
           'Closed Constraints file has wrong number of inputs (evalFunction)');
end
if ~feasible
   iterate.c = [];
   iterate.h = Inf;
   iterate.f = Inf;
   if grad
      iterate.gradc = [];
      iterate.gradh = Inf;
      iterate.gradf = Inf;
   end
end

%*******************************************************************************
% evalFunction:  Evaluates objective and constraint functions at a point.
% ------------------------------------------------------------------------------
% Called by: evalPointSet
% Calls: < User F >
% VARIABLES:
%  iterate  = current iterate
%    .x     =   coordinates of iterate
%    .p     =   values of categorical variables
%    .f     =   objective function value of iterate
%    .c     =   constraint function values of iterate
%    .h     =   constraint violation function value of iterate
%    .gradf =   f'(x) at current iterate x
%    .gradc =   c'(x) at current iterate x
%    .gradh =   h'(x) at current iterate x
%  f        = name of optimization problem functions file
%  ftype    = type of functions file (F=FORTRAN, C=C/C++, M=MATLAB)
%  grad     = flag for computing gradients, if available
%  feasible = logical indicating if closed constraints are satisfied
%  pchar    = categorical variables stored as character strings
%  pint     = categorical variables stored as integers
%  cflag    = flag for sorting out char and int categorical variables
%  nc       = number of nonlinear constraints
%  badValue = the value used when f(x) is complex (Inf)
%  cplus    = vector of nonlinear constraint violations
%*******************************************************************************
function iterate = evalFunction(f,ftype,iterate,grad)

% Compute surrogate h(x) (and h'(x) if available), if appropriate
if isfield(iterate,'sc')
   cplus = (iterate.sc > 0).*iterate.sc;
   iterate.sh = norm(cplus)^2;
   if grad
      iterate.sgradh = 2*iterate.sgradc*cplus';
   end
end

% Evaluate functions file
switch upper(ftype)

% Process a compiled Fortran functions file
case {'F','C'}
   if isempty(iterate.p)
      pchar = 0;
      pint  = 0;
   else
      cflag = zeros(length(iterate.p),1);
      for k = 1:length(iterate.p)
         cflag(k) = ischar(iterate.p{k});
      end
      pchar = [iterate.p{find( cflag)}];  %#ok
      pint  = [iterate.p{find(~cflag)}];  %#ok
   end

   [nc,iterate.f,c,gradf,gradc] = feval(f,iterate.x,pint,pchar);
   iterate.c = c(1:nc);
   if grad
      iterate.gradf = gradf;
      iterate.gradc = gradc(:,1:nc);
   end

% Process a Matlab functions file
case {'M'}
   switch nargin(f) + 2*grad
   case {1}
      if nargout(f) == 1
         iterate.f = feval(f,iterate.x);
         iterate.c = [];
      else
         [iterate.f,iterate.c] = feval(f,iterate.x);
      end
   case {2}
      if nargout(f) == 1
         iterate.f = feval(f,iterate.x,iterate.p);
         iterate.c = [];
      else
         [iterate.f,iterate.c] = feval(f,iterate.x,iterate.p);
      end
   case {3}
      [iterate.f,iterate.c,iterate.gradf,iterate.gradc] = feval(f,iterate.x);
   case {4}
      [iterate.f,iterate.c,iterate.gradf,iterate.gradc] = feval(f,iterate.x, ...
                                                                  iterate.p);
   otherwise
      error('mads:function:badVariables', ...
            'Error in EvalFunction logical input variables (evalFunction).');
   end
otherwise
   error('mads:function:badType','Invalid Function Type (evalFunction).');
end

% Prevent complex-valued functions
badValue = realmax;
iterate.f(isempty(iterate.f)) = badValue;
iterate.f(~isreal(iterate.f)) = badValue;
iterate.c(~isreal(iterate.c)) = badValue;
if grad
   iterate.gradf(~isreal(iterate.gradf)) = badValue;
   iterate.gradc(~isreal(iterate.gradc)) = badValue;
end

% Compute h(x) (and h'(x) if available)
cplus = (iterate.c > 0).*iterate.c;
iterate.h = norm(cplus)^2;
if grad
   iterate.gradh = 2*iterate.gradc*cplus';
end

% Append to the Cache any user-specified parameter
if isappdata(0,'PARAM')
   Param = getappdata(0,'PARAM');
   if isfield(Param,'param')
      iterate.param = Param.param;
   else
      iterate.param = 0;
   end
end
return

%*******************************************************************************
% updateCache:  Update the Cache with the current iterate.
% ------------------------------------------------------------------------------
% Called by: EvalPointSet
% VARIABLES:
%  Cache      = structure containing data on all past iterates
%    .iterate =   the previously computed iterates
%    .size    =   current number of iterates in the Cache
%    .tol     =   tolerance for determining if iterate was in Cache
%    .xnorm   =   vector of 1-norms of previously computed iterates
%    .pID     =   vector of categorical variable ID numbers
%    .bfp     =   vector of best feasible points
%    .lip     =   vector of least infeasible points
%  Problem    = structure containing problem file names and parameters
%    .isMVP   =   logical indicating if problem is an MVP
%    .maxNc   =   maximum number of constraints to be added
%    .maxNx   =   maximum number of continuous variables of any iterates
%    .maxNp   =   maximum number of categorical variables of any iterates
%  iterate    = current iterate
%    .x       =   coordinates of iterate
%    .p       =   values of categorical variables
%    .n       =   number of continuous variables
%    .f       =   objective function value of iterate
%    .c       =   constraint function values of iterate
%    .h       =   constraint violation function value of iterate
%    .gradf   =   f'(x) at current iterate x
%    .gradc   =   c'(x) at current iterate x
%    .gradh   =   h'(x) at current iterate x
%  xnorm      = 1-norm of the curent iterate
%  pID        = categorical variable value ID of the current iterate
%  nFunc      = number of function evaluations at point (=1, if not stochastic)
%  sdev       = standard deviation at point (=0, if not stochastic)
%  n          = number of iterates for which memory is to be allocated
%  ind        = index variable for additional chunks of memory
%  p          = temporary storage
%*******************************************************************************
function Cache = updateCache(Cache,Problem,iterate,xnorm,pID,nFunc,sdev)

% Allocate extra chunks of memory if more is required
if Cache.size >= max(1024,length(Cache.iterate))
   n = 1024;
   ind = length(Cache.iterate) + (1:n);
   if Problem.isMVP
      [p{1:Problem.maxNp}] = deal('            ');
   else
      p = {};   
   end
   [Cache.iterate(ind).x]     = deal(zeros(Problem.maxNx,1));
   [Cache.iterate(ind).p]     = deal(p);
   [Cache.iterate(ind).n]     = deal(0);
   [Cache.iterate(ind).type]  = deal(' ');
   if isfield(iterate,'sf')
      [Cache.iterate(ind).sf] = deal(0);
   end
   if isfield(iterate,'sc')
      [Cache.iterate(ind).sc] = deal(zeros(Problem.maxNc,1));
   end
   if isfield(iterate,'sh')
      [Cache.iterate(ind).sh] = deal(0);
   end
   if isfield(iterate,'sgradf')
      [Cache.iterate(ind).sgradf] = deal(zeros(Problem.maxNx,1));
   end
   if isfield(iterate,'sgradc')
      [Cache.iterate(ind).sgradc] = deal(zeros(Problem.maxNc,Problem.maxNx));
   end
   if isfield(iterate,'sgradh')
      [Cache.iterate(ind).sgradh] = deal(zeros(Problem.maxNx,1));
   end
   [Cache.iterate(ind).f]     = deal(0);
   [Cache.iterate(ind).c]     = deal(zeros(Problem.maxNc,1));
   [Cache.iterate(ind).h]     = deal(0);
   [Cache.iterate(ind).param] = deal(zeros(1,length(iterate.param)));
   if isfield(iterate,'gradf')
      [Cache.iterate(ind).gradf] = deal(zeros(Problem.maxNx,1));
      [Cache.iterate(ind).gradc] = deal(zeros(Problem.maxNc,Problem.maxNx));
      [Cache.iterate(ind).gradh] = deal(zeros(Problem.maxNx,1));
   end
   Cache.xnorm(ind)           = deal(0);
   Cache.pID(ind)             = deal(0);
   Cache.bfp(ind)             = deal(0);
   Cache.lip(:,ind)           = deal(0);
   Cache.nFunc(ind)           = deal(0);
   Cache.sdev(ind)            = deal(0);
   Cache.isSurrPoint(ind)     = deal(0);
end

% Update Cache (.isSurrPoint is set to zero and updated later)
Cache.size                    = Cache.size + 1;
Cache.iterate(Cache.size)     = orderfields(iterate,Cache.iterate);
Cache.xnorm(Cache.size)       = xnorm;
Cache.pID(Cache.size)         = pID;
Cache.nFunc(Cache.size)       = nFunc;
Cache.sdev(Cache.size)        = sdev;
Cache.isSurrPoint(Cache.size) = 0;

% Update the BFP and LIP function values in the Cache (used in plotting)
feasible = (iterate.h <= Cache.Filter.hmin);
if feasible
   if Cache.size == 1
      Cache.bfp(Cache.size)   = iterate.f;
      Cache.lip(:,Cache.size) = [Inf; Inf];
   else
      Cache.bfp(Cache.size)   = min(Cache.bfp(Cache.size-1),iterate.f);
      Cache.lip(:,Cache.size) = Cache.lip(:,Cache.size-1);
   end
else
   if Cache.size == 1
      Cache.bfp(Cache.size)   = Inf;
      Cache.lip(:,Cache.size) = [iterate.f; iterate.h];
   else
      Cache.bfp(Cache.size) = Cache.bfp(Cache.size-1);
      if iterate.h < Cache.lip(Cache.size-1)
         Cache.lip(:,Cache.size) = [iterate.f; iterate.h];
      else
         Cache.lip(:,Cache.size) = Cache.lip(:,Cache.size-1);
      end
   end
end

return

%*******************************************************************************
% FUNCTIONS FOR UPDATING AND PLOTTING THE FILTER
%*******************************************************************************
%*******************************************************************************
% updateFilter:  Update the filter and solutions vectors, given in iterate.
% ------------------------------------------------------------------------------
% Called by: evalPointSet, mvpPoll
% Calls:     dominates, plotFilter
% VARIABLES:
%  unfiltered  = logical indicating that the iterate is unfiltered
%  Filter      = structure containing filter
%    .hmin     =   minimum allowed constraint violation of filter point
%    .hmax     =   maximum allowed constraint violation of filter point
%    .strict   =   flag indicating that hmax is strictly enforced
%    .plot     =   flag for displaying a real-time plot of the filter 
%    .feasible =   sorted indices of feasible Cache iterates
%    .F        =   sorted indices of Cache iterates in the filter
%  iterate     = current iterate
%    .x        =   coordinates of iterate
%    .p        =   values of categorical variables
%    .n        =   dimension of the continuous variables
%    .f        =   objective function value of iterate
%    .c        =   constraint function values of iterate
%    .h        =   constraint violation function value of iterate
%    .gradf    =   f'(x) at current iterate x
%    .gradc    =   c'(x) at current iterate x
%    .gradh    =   h'(x) at current iterate x
%  Cache       = collection of all previously evaluated iterates
%    .iterate  =   vector of iterates
%    .size     =   number of iterates
%  maxSolPts   = maximum number of feasible points to keep in filter
%  ind         = indices used in sorting
%  dominated   = flags indicating which filter points iterate dominates
%  infeasible  = flags indicating which filter points exceed hmax
%*******************************************************************************
function [unfiltered,Filter] = updateFilter(iterate,Filter,Cache)

maxSolPts = 10;
feasible  = iterate.h <= Filter.hmin;

% FEASIBLE CASE: Update sorted list of 10 best feasible solutions
if feasible
   fvalues = [Cache.iterate(Filter.feasible).f];
   ind     = find(fvalues > iterate.f,1);
   if isempty(ind)
      if isempty(Filter.feasible), Filter.feasible(1) = Cache.size; end
   else
      Filter.feasible = [Filter.feasible(1:ind-1); Cache.size; ...
                         Filter.feasible(ind:end)];
      Filter.feasible(maxSolPts:end) = [];
   end
   Filter.BFP = Filter.feasible(1);
   unfiltered = (Filter.BFP == Cache.size);
   return
end

% INFEASIBLE CASE

% Initialize hmax and test to see if point is unfiltered
if ~Filter.strict && isempty(Filter.F)
   hmax = 10*iterate.h;
elseif ~Filter.strict && Cache.iterate(Filter.F(1)).h >= Filter.hmax
   hmax = 10*Cache.iterate(Filter.F(1)).h;
else
   hmax = Filter.hmax;
end
unfiltered = iterate.h<hmax && ~any(dominates(Cache.iterate(Filter.F),iterate));

% Delete newly dominated points, add the unfiltered point, and replot
if unfiltered
   dominated  = dominates(iterate,Cache.iterate(Filter.F));
   infeasible = [Cache.iterate(Filter.F).h] >= hmax;
   Filter.F(dominated | infeasible) = [];
   hvalues = [Cache.iterate(Filter.F).h];
   ind     = find(hvalues > iterate.h,1);
   if isempty(ind)
      Filter.F(end+1) = Cache.size;
   else
      Filter.F = [Filter.F(1:ind-1), Cache.size, Filter.F(ind:end)];
   end
   Filter.LIP = Filter.F(1);
   Filter.IIP = Filter.F(end);
   if Filter.plot, plotFilter(Filter,hmax,Cache.iterate,maxSolPts); end
end
return

%*******************************************************************************
% dominates:  Determines if x dominates y with respect to .f and .h values.
% ------------------------------------------------------------------------------
% Called by: updateFilter
% VARIABLES:
%  d    = logicals indicating whether or not x dominates y
%  x,y  = two iterates to be compared
%    .f =   objective function value
%    .h =   aggregated constraint violation function value
%*******************************************************************************
function d = dominates(x,y)
d = ~([y.f] < [x.f] | [y.h] < [x.h]);
return

%*******************************************************************************
% plotFilter:  Plot the filter.
% ------------------------------------------------------------------------------
% Called by: updateFilter, createCache
% VARIABLES:
%  Filter        = structure containing the filter
%    .F          = sorted indices of Cache iterates in the filter
%    .feasible   = sorted indices of feasible Cache iterates
%    .hmax       = maximum allowed constraint violation of filter point
%    .plothandle = plot handle
%  hmax          = the filter hmax for the plot only
%  Cache         = collection of all previously evaluated iterates
%    .f          =   objective function value
%    .h          =   constraint violation function value
%  nPoints       = number of filter points to plot
%  nFilter       = number of iterates in the filter
%  nSolutions    = number of feasible solutions
%  nPlot         = number of feasible solutions to plot
%  feasible_f    = vector of f-values of feasible points to be plotted
%  h             = vector of h-values of filter points to be plotted
%  f             = vector of f-values of filter points to be plotted
%  pctAxis       = percentage of extra space between axes and data
%  ax            = plot axis parameters
%*******************************************************************************
function plotFilter(Filter,hmax,Cache,nPoints)

if isempty(Filter.F), return; end
nSolutions = length(Filter.feasible);
nPlot      = min(nPoints,nSolutions);

% Create feasible solutions vector and filter vectors for plotting
feasible_f = [Cache(Filter.feasible(:)).f];
[h,f] = stairs([Cache(Filter.F).h],[Cache(Filter.F).f]);
h = [h(1); h];
if nPlot
   f = [max(Cache(Filter.feasible(nPlot)).f,Cache(Filter.F(1)).f); f];
else
   f = [Cache(Filter.F(1)).f; f];
end
h(end+1) = hmax;
f(end+1) = f(end);

% Create filter plot
set([Filter.plothandle; get(Filter.plothandle,'Children')],'Visible','on');
if nPlot
   plot(Filter.plothandle,h,f,'k-',zeros(nPlot,1),feasible_f(1:nPlot),'b*');
else
   plot(Filter.plothandle,h,f,'k-',[],[],'b*');
end
pctAxis = 0.02;
ax = axis;
ax(1:2) = ax(1:2) + pctAxis*(ax(2)-ax(1))*[-1,1];
ax(3:4) = ax(3:4) + pctAxis*(ax(4)-ax(3))*[-1,1];
line([h(end),h(end)],[f(end) ax(3)]);
axis(ax);
title(Filter.plothandle, 'Filter','FontWeight', 'bold','FontSize',11);
xlabel(Filter.plothandle,'Constraint Violation, h(x)', 'FontSize',12);
ylabel(Filter.plothandle,'Objective Function, f(x)  ', 'FontSize',12);
drawnow;
return


%*******************************************************************************
% FUNCTIONS FOR TERMINATION
%*******************************************************************************
%*******************************************************************************
% processOutput:  Process output and delete temporary files and variables.
% ------------------------------------------------------------------------------
% Called by: mads
% Calls:     getPollCenter, plotHistory
% VARIABLES:
%  BestF           = best feasible point found
%  BestI           = least infeasible point found
%  Problem         = structure containing optimization problem data
%    .File.P       =   name of parameter file
%    .nameCache    =   name of the base workspace Cache variable
%  Options         = structure of MADS parameters
%    .plotHistory  =   flag for plotting f-value vs function evals
%    .hplothandle  =   handle of the history plot
%    .plotColor    =   string indicating color of history plot line
%    .saveHistory  =   flag for saving history to file
%  RunData         = structure of run performance data
%    .iCache       =   indices of all Cache points that were poll centers
%  Cache           = structure containing previously evaluated iterates
%    .Filter       =   structure containing filter data
%      .hmin       =     minimum allowable h value to be in Filter.F
%  nFunc           = vector of cumulative numbers of function evaluations
%  ICache          = partial Cache of only poll centers
%    .nFunc        =  cumulative number of function evaluations
%    .f            =  objective function value
%    .x            =  continuous variable values
%  junk            = temporary unused stroage
%  fType           = type of file to be saved (.txt or .mat)
%*******************************************************************************
function [BestF,BestI] = processOutput(Problem,Options,RunData)

% Plot performance history
Cache = getappdata(0,Problem.nameCache);
if Options.plotHistory
   plotHistory(Options.hplothandle,Options.plotColor,Cache);
end

% Save to file and plotting
if Options.saveHistory
   nFunc  = cumsum(Cache.nFunc);
   nFunc  = nFunc(RunData.iCache);
   ICache = Cache.iterate(RunData.iCache);
   for k = 1:length(ICache)
      ICache(k).nFunc = nFunc(k);
   end
   [~,fType] = strtok(Problem.File.H,'.');
   if strcmp(fType,'mat')
      save(Problem.File.H, 'ICache','-V6');
   else
      fid = fopen(Problem.File.H,'at');
      for k = 1:length(ICache)
         fprintf(fid,'%-2s %6d %s %s %s \n', ICache(k).type, ...
                 ICache(k).nFunc, num2str(ICache(k).f,15), ...
                 num2str(ICache(k).x',12), int2str(cell2mat(ICache(k).p)'));
      end
      fclose(fid);
      rehash;
   end
end

% Store Best Solutions
Best = getPollCenter([0,1],Cache.Filter,Cache,Options.hRho);
if length(Best) == 1
   [BestF,BestI] = deal(Best);  
else
   [BestF,BestI] = deal(Best(1),Best(2));
end
if BestF.h > Cache.Filter.hmin, BestF = []; end
if BestI.h < Cache.Filter.hmin || BestI.h == 0, BestI = []; end
return

%*******************************************************************************
% closeWorkspace:  Shuts down base workspace appdata and deletes temp files.
% ------------------------------------------------------------------------------
% Called by: mads
% VARIABLES:
%  Cache           = structure containing previously evaluated iterates
%  Problem         = structure containing optimization problem data
%    .File.F       =   name of functions file
%    .nameCache    =   name of the base workspace Cache variable
%  RunData         = structure containing MADS run statistics
%  surrName        = string used in constructing surrogate problem filename
%*******************************************************************************
function Cache = closeWorkspace(Problem)

% Store Cache and delete all base workspace appdata
Cache = getappdata(0,Problem.nameCache);
if isappdata(0,Problem.nameCache), rmappdata(0,Problem.nameCache); end
if isappdata(0,'ALTCONE'),         rmappdata(0,'ALTCONE');         end
if isappdata(0,'SUR') && strcmp(Problem.nameCache,'CACHE')
   rmappdata(0,'SUR');
end
if isappdata(0,'PS') && strcmp(Problem.nameCache,'CACHE')
   rmappdata(0,'PS');
end

% Delete function files used in search and surrogate optimization
surrName = [Problem.File.F, '_Sur'];
lastwarn('');
% if exist('SurMVPOmega',      'file'), delete('SurMVPOmega.m');       end
if exist([surrName,'_Cache'],'file'), delete([surrName,'_Cache.m']); end
if ~isempty(lastwarn)
   warning(['Search or surrogate file found but not deleted.      \n', ...
            'File may have been moved from its original location. \n', ...
            'Please delete the file manually.'],[]);
end
if strcmpi(Problem.typeProblem,'TRUTH')
   diary off
end
return

%*******************************************************************************
% plotHistory:  Plot optimal value vs. number of function evaluations.
% ------------------------------------------------------------------------------
% Called by: update, ProcessOutput
% VARIABLES:
%  handle  = handle of the plot axes
%  color   = string code indicating color of plot curve
%  symbol  = string code indicating linestyle/symbol of plot curve
%  Cache   = Cache of previously computed iterates
%    .bfp  =   vector of best feasible f-values
%    .lip  =   vector of least infeasible f- and h-values
%    .size =   number of iterates in the Cache
%*******************************************************************************
function plotHistory(handle,color,Cache)

% Set up history plot axes handle
if isempty(handle) || ~ishandle(handle)
   figure; handle = gca;
end

% Construct f plot using outside axes handle, if avaliable
set([handle; get(handle,'Children')],'Visible','on');
if Cache.size <= 1
   symbol = '*';
else
   symbol = '-';
end
if any(isfinite(Cache.bfp))
   plot(handle,cumsum(Cache.nFunc(1:Cache.size),2), ...
               Cache.bfp(1:Cache.size)',  [color,symbol],'LineWidth',2);
else
   plot(handle,cumsum(Cache.nFunc(1:Cache.size),2),...
               Cache.lip(1,1:Cache.size)',[color,symbol],'LineWidth',2);
end
title(handle, 'Performance History', 'FontWeight','bold', 'FontSize',11);
xlabel(handle,'Number of Function Evaluations',           'FontSize',12);
ylabel(handle,'Objective Function Value  ',               'FontSize',12);
drawnow;
return

%*******************************************************************************
% terminate:  Tests satisfaction of termination criteria
% ------------------------------------------------------------------------------
% Called by: mads
% VARIABLES:
%  term         = flag indicating if program should be terminated
%  R            = structure of run data
%    .delta     =   current poll size
%    .nIter     =   current iteration number
%    .nFunc     =   current number of function evaluations
%    .time      =   current cumulative CPU time
%    .nFails    =   current number of consecutive poll failures
%    .stopRun   =   current value of the Stop Run button
%  T            = structure of termination criteria
%    .delta     =   termination poll size
%    .nIter     =   maximum iteration number
%    .nFunc     =   maximum number of function evaluations
%    .time      =   maximum cumulative CPU time
%    .nFails    =   maximum number of consecutive poll failures
%*******************************************************************************
function term = terminate(T,R)

% Set up successful termination conditions
term = R.delta <= T.delta;

% Add in unsuccessful termination conditions
term = term || R.nIter >= T.nIter || R.nFunc  >= T.nFunc  || ...
               R.time  >= T.time  || R.nFails >= T.nFails || R.stopRun;
return

%*******************************************************************************
% UTILITIES
%*******************************************************************************
%*******************************************************************************
% rethrow2: Alternate error rethrowing function
% ------------------------------------------------------------------------------
% Called by: runMADS
% VARIABLES:
%   err        = structure of error data
%     .message = error message
%     .stack   = error stack
%       .file  = full filename(s) of file(s) in which error occurred
%       .name  = names(s) of function(s) in which error occurred
%       .line  = line number(s) where error occurred
%   fname      = filename(s) without path
%   name       = file/function string where error occurred
%*******************************************************************************
function rethrow2(err)

disp(['??? ', err.message]);
disp(' ');
for k = 1:length(err.stack)
   [~,fname] = fileparts(err.stack(k).file);
   if strcmp(fname,err.stack(k).name)
      name = fname;
   else
      name = [fname,'>',err.stack(k).name];
   end
   disp(['Error in ==> ',name,' at ',int2str(err.stack(k).line)]);
end
disp(' ');
if (1), rethrow(err); end
return  %#ok

%*******************************************************************************
% runLP:  Runs a linear programming (LP) problem
% ------------------------------------------------------------------------------
% Called by: makeFeasible, gradDirections, removeRedundancy, updateOmega
% Calls:     linprog (Optim_TB), optimset(Optim_TB)
% VARIABLES:
%  x    = solution to LP problem
%  fval = objective function value at the solution to the LP problem
%  err  = error flag for LP
%  c    = vector objective function coefficients
%  A    = matrix of inequality constraint coefficients
%  b    = vector of inequality constraint right-hand sides
%  Aeq  = matrix of equality constraint coefficients
%  beq  = vector of equality constraint right-hand sides
%  LU   = vector of lower bounds
%  UB   = vector of upper bounds
%  x0   = initial point
%*******************************************************************************
function [x,fval,err] = runLP(c,A,b,Aeq,beq,LB,UB,x0)

% Test to see if Optimization Toolbox license is present
if license('test','Optimization_Toolbox')

   % Handle missing arguments
   if nargin < 8, x0  = []; end
   if nargin < 7, UB  = []; end
   if nargin < 6, LB  = []; end
   if nargin < 5, beq = []; end
   if nargin < 4, Aeq = []; end

   % Solve LP problem
   warning('off','all');
   [x,fval,err] = linprog(c,A,b,Aeq,beq,LB,UB,x0,optimset('Display','off'));
   warning('on','all');
   
% Default values if license not present
else
   x = []; fval = []; err = -9;
end
return

%*******************************************************************************
% bisect:  Bisection 1-D root-finding algorithm with rounding option.
% ------------------------------------------------------------------------------
% Called by: adjHalton, daceSurrogate
% VARIABLES:
%  x    = final solution (and midpoint between shrinking subintervals)
%  fx   = function value at final solution
%  f    = function for which a root will be found
%  a    = lower bound on interval containing root
%  b    = upper bound on interval containing root
%  tol  = error tolerance on solution
%  flag = flag for rounding final solution to left or right {-1,1} endpoint
%*******************************************************************************
function [x,fx] = bisect(f,a,b,tol,flag)

fa = f(a);
fb = f(b);
if fa*fb > 0
   error('MADS:bisect','No root detected in interval');
end
while (b - a)/2 > tol
   x = (a+b)/2;
%   x = a - fa*(b-a)/(fb-fa);
   fx = f(x);
   if fx == 0.0, return, end
   if fa*fx < 0
      b  = x;
   else
      a  = x;
      fa = fx;
   end
end
x = x + flag*(b-a)/2;
return

%*******************************************************************************
% cat2cont:  Converts a mixed variable iterate to a continuous iterate.
% ------------------------------------------------------------------------------
% Called by: search, getDataSites, optimizeSurrogate, inOmega
% VARIABLES:
%  iterate  = the iterate to be converted
%  plist    = cell array of lists of possible catergorical variable values
%*******************************************************************************
function [iterate] = cat2cont(iterate,plist)

for k = 1:length(iterate.p)
   if ischar(iterate.p{k})
      iterate.x = [iterate.x; find(strcmp(plist{k},iterate.p{k}))];
   else
      iterate.x = [iterate.x; iterate.p{k}];
   end
end
iterate.n = length(iterate.x);
return

%*******************************************************************************
% dealIterate:  Assigns columns of a matrix to a vector of iterates.
% ------------------------------------------------------------------------------
% Called by: search, poll, optimizeSurrogate
% VARIABLES:
%  S       = vector of nS iterates
%  iterate = sample iterate with n continuous variables
%  nS      = number of iterates
%  X       = n x nS matrix
%  plist   = cell array of allowable categorical variable values
%  np      = number of categorical variables
%*******************************************************************************
function S = dealIterate(iterate,nS,X,varargin)

if nargin > 3
   plist = varargin{1};
else
   plist = {};
end
n  = length(iterate.x);
np = length(plist);
if all(size(X) == [n+np,nS])
   [S(1:nS).x] = deal(iterate.x);
   [S(1:nS).p] = deal(iterate.p);
   for k = 1:nS
      S(k).x = X(1:n,k);
      for j = 1:np
         S(k).p{j} = plist{j}{X(n+j,k)};
      end
   end
else
   error('dealIterate:Input','Invalid matrix dimensions');
end
return

%*******************************************************************************
% myeval: Evaluate variable of unknown type.
% ------------------------------------------------------------------------------
% Called by: cmaes_mod
% VARIABLES:
%  s     = variable of unknown type
%  res   = result after evaluating s
%  fname = cell array of fieldnames of s if s is a structure
%*******************************************************************************
function res = myeval(s)

if isstruct(s)
   fname = fieldnames(s);
   for k = 1:length(fname)
      if ischar(s.(fname{k}))
         res.(fname{k}) = evalin('caller', s.(fname{k}));
      else
         res.(fname{k}) = s.(fname{k});
      end      
   end
elseif iscell(s)
   for k = 1:length(s)
      if ischar(s{k})
         res(k) = evalin('caller', s{k}); %#ok
      else
         res(k) = s{k}; %#ok
      end
   end
elseif ischar(s)
   res = evalin('caller', s);
else
   res = s;
end
return

%*******************************************************************************
% default_N:  Default discrete neighbor file for MVP problems
% ------------------------------------------------------------------------------
% Called by: mvpPoll
% VARIABLES:
%  N        = set of discrete neighbors
%  iterate  = current iterate
%  plist    = list of allowed values that categorical variables may take on
%  neighbor = constructed neighboring point of the current iterate
%  ind      = index into plist where a categorical variable value matches  
%*******************************************************************************
function N = default_N(Problem,iterate,plist,delta)  %#ok

N = [];
for k = 1:length(iterate.p)
   if ischar(iterate.p{k})
      ind    = find(strcmp(iterate.p{k},plist{k}));
      nplist = length(plist{k});
   else
      ind = find([plist{k}{:}] == iterate.p{k});
      nplist = length(plist{k}{:});
   end
   ind = ind(1);
   if ind > 1
      neighbor.x = iterate.x;
      neighbor.p = iterate.p;
      if ischar(iterate.p{k})
         neighbor.p{k} = plist{k}{ind-1};
      else
         temp = [plist{k}{:}];
         neighbor.p{k} = temp(ind-1);
      end
      N = [N, neighbor];  %#ok
   end
   if ind < nplist
      neighbor.x = iterate.x;
      neighbor.p = iterate.p;
      if ischar(iterate.p{k})
         neighbor.p{k} = plist{k}{ind+1};
      else
         temp = [plist{k}{:}];
         neighbor.p{k} = temp(ind+1);
      end
      N = [N, neighbor];  %#ok
   end
end
return

%*******************************************************************************
% dec2bigbase:  Convert decimal integer to base B vector
% ------------------------------------------------------------------------------
% Called by: adjustedHalton
% VARIABLES:
%*******************************************************************************

%DEC2BIGBASE Convert decimal integer to base B vector.
%   DEC2BIGBASE(D,B) returns the representation of D as a vector of
%   digits in base B.  D must be a non-negative integer smaller than 2^52
%   and B must be an integer greater than 1.  
%
%   DEC2BIGBASE(D,B,N) produces a representation with at least N digits.
%
%   Examples
%       dec2bigbase(23,3) returns [2 1 2] 
%       dec2bigbase(23,3,5) returns [0 0 2 1 2]
%
%   See also DEC2BASE, BASE2DEC, DEC2HEX, DEC2BIN.

%   written by Douglas M. Schwarz
%   Eastman Kodak Company (on leave until 4 Jan 1999)
%   schwarz@kodak.com, schwarz@servtech.com
%   1 October 1998
function s = dec2bigbase(d,base,n)

 narginchk(2,3);

 d = d(:);

 base = floor(base);
 if base < 2, error('B must be greater than 1.'); end
 if base == 2
    [~,nreq] = log2(max(d));
 else
    nreq = ceil(log2(max(d) + 1)/log2(base)); 
 end

 if nargin == 3
    nreq = max(nreq,1);
    n    = max(n,nreq);
    last = n - nreq + 1;
 else
    n    = max(nreq,1);
    last = 1;
 end

 s(:,n) = rem(d,base);
 while n ~= last
    n = n - 1;
    d = floor(d/base);
    s(:,n) = rem(d,base);
 end
return


%===============================================================================
% FUNCTIONS ADAPTED FROM THE CMAES MATLAB CODE BY NIKOLAUS HANSEN
%===============================================================================
%===============================================================================
% cmaes_mod: Evolutionary strategy with covariance matrix adaptation for solving
%            optimization problems.
%   Adapted for NOMADm from the CMA-ES MATLAB package by Nikolaus Hansen
% ------------------------------------------------------------------------------
% VARIABLES
%  xmin                = minimal point found from last iteration
%  fmin                = function value at xmin
%  nFunc               = number of function evaluations performed
%  out                 = structure of results from all previous runs
%  f                   = string name or function handle of objective function
%  x0                  = initial point
%  sigma0              = initial standard deviations for each variable
%  Options             = structure of user-defined options
%    .Stop.fitness     =   threshold for stopping if f(x) falls below it
%    .Stop.maxFunEval  =   maximum number of function evaluations
%    .Stop.maxIter     =   maximum number of iterations
%    .Stop.tolX        =   smallest change in x allowed
%    .Stop.tolUpX      =   largest change in x allowed
%    .Stop.tolFun      =   smallest change in F(x) allowed
%    .debug            =   numeric code indicating level of debugging
%    .diffMaxChange    =   maximal variable changes (vector/scalar)
%    .diffMinChange    =   minimal variable changes (vector/scalar)
%    .lBounds          =   lower bounds (vector/scalar)
%    .uBounds          =   upper bounds (vector/scalar)
%    .evalInitialX     =   flag indicating that given initial is to be evaluated
%    .incPopSize       =   population size multiplier before restart
%    .popSize          =   population size, lambda
%    .nParents         =   popsize = lambda
%    .nRestarts        =   number of restarts
%    .recombWeights    = type of recombination weighting
%    .science          =   flag for additional (minor) problem capturing (=off)
%    .seed             =   random number generator seed
%  varargin            = additional arguments passed to objective function 
%  cc                  = time constant for cumulation for covariance matrix
%  cs                  = time constant for cumulation for step size control
%  ccov                = learning rate for covariance matrix
%  mucov               = size of mu used for calculating learning rate ccov
%  damps               = damping for controlling sigma (set ~1, limits increase)
%  sigma               = overall standard deviation
%  B                   = coordinate axes
%  C                   = covariance matrix
%  D                   = diagonal scaling matrix
%  BD                  = B*D
%  pc                  = evolution path for C
%  ps                  = evolution path for sigma
%  fitness             = structure of past fitness values
%    .hist             =   vector of fitness values
%    .histsel          =   vector of adjusted fitness values
%  bnd                 = structure of variable bound information
%    .exists           =   flag indicating that at least one bound exists
%    .weights          =   weights for penalizing infeasibility
%  stopflag            = cell array of flags set when termination occurred
%  bestever            = structure containing minimal point across all runs
%    .x                =   variable values
%    .f                =   function value
%    .nFunc            =   number of function evaluations
%-------------------------------------------------------------------------------
% OTHER COMMENTS:
% 1.  For hard nonlinear constraints, f must return a NaN value at infeaible 
%     points, in which case, the function evaluation will not be counted.
%
% 2.  The initial x0 point can be a vector, matrix, or string.  If X0 is a 
%     matrix, then mean(X0,2) is taken as the initial point.  If X0 is a string,
%     such as '2*rand(10,1)-1', the initial point is generated by evaluating the
%     string.
%
% 3.  The input variable sigma0 can be a scalar, column vector or string to be
%     evaluated.  It determines the initial standard deviations for each
%     variable.  An appropriate value for each variable would be 1/3 of the
%     variable's range.  For example, for a range of [0,6]^10, sigma0 = 2, in
%     which case, the call to cmaes_mod would be cmaes('myfun',3*rand(10,1),2).
%     The default if x0 is a matrix is sqrt(var(x0')').
%
% 4.  The original default values for Options have been moved to the calling
%     function, so that unspecified values are assigned their default values.
%     In many cases, a string value is evaluated as an expression.
%
% 5.  Possible values for stopflag are 'fitness', 'tolx', 'tolupx', 'tolfun',
%     'maxfunevals', 'maxiter', 'manual', 'bug'.  To make use of these for
%     'tolx', for example, use the expressions, any(strcmp(STOPFLAG, 'tolx'))
%     or findstr(strcat(STOPFLAG,'tolx')).
%
% 6.  The structure out contains results of the runs.  Evolution of the mean 
%     value is contained in out.hist.mean and the overall best solution is
%     found in out.solutions.bestever.  Each of these has subfields, .x, .f. and
%     .evals, that represent the point, it's function value, and the number of
%     function evaluations performed.
%
% 7.  A manual stop is performed if the first two non-blank sequences in any
%     line of the file, ABORT.txt,  are 'stop'.  This is useful for problems
%     with computationally expensive functions.  
%
% 8.  A primary tunable parameter is the popluation size, Options.popSize 
%     (lambda).  Increasing this value, which by default also increases
%     Options.nParents, improves the global nature of the solution in
%     exchange for computational effort (i.e., speed).  As a general rule, speed
%     decreases at most linearly as population size is increased.  The authors
%     recommend a strategy of increasing lambda by a factor of 3, starting with
%     initial values between 5 and 20.  They assert that even sizes as large as 
%     1000+100*N can be reasonable.
%
% 10. Options.nRestarts and Options.incPopSize can be used to automate a
%     multistart procedure, where the PopSize is increased by a factor of
%     IncPopSize (default = 2) before each restart.  For each restart, x0 (given
%     as a string) is reevaluated.
%
% 12. The adaptation of the covariance matrix (e.g., by the CMA) is equivalent
%     to a general linear transformation of the problem coding.  Nevertheless, 
%     all problem-specific knowledge about the best linear transformation should
%     be exploited before starting the search.  That is, an appropriate a priori 
%     transformation should be applied to the problem, ideally so that the
%     identity matrix would make the best choice as the initial covariance
%     matrix.
%
% 13. For lower and upper bounds, use Options.lBound and Options.uBound.  They
%     can be given as scalars or column vectors.
%===============================================================================
function [xmin,fmin,nFunc,out] = cmaes_mod(f,x0,sigma0,Options,varargin)

% Process input parameters
if nargin < 1 || isempty(f) || ~(isa(f,'function_handle') || ischar(f))
   error('cmaes:f','First argument F must be nonempty string');
end
if nargin < 2 || isempty(x0)
   error('cmaes:x0','Second argument X0 must be nonempty vector');
end
if nargin < 4 || isempty(Options)
   error('cmaes:options','Options must be defined');
end
if nargin < 3 || isempty(sigma0)
   if size(myeval(x0),2) > 1
      sigma0 = std(x0,0,2);
      if any(~sigma0)
         error('cmaes:sigma', ...
               'Initial search volume is zero. Choose SIGMA or X0 appropriate');
      end
   else
      sigma0 = [];
   end
end

% Initialize data storage
out.Cache.x = [];
out.Cache.f = [];

% Begin main loop
stoplabels = {'bug','fitness','maxfunevals','maxiter','tolx','tolupx',...
              'tolfun','manual'};
nFunc    = 0;
for irun = 1:myeval(Options.nRestarts) + 1

   % Assign settings from input parameters and options for myeval
   xmean   = mean(myeval(x0),2);
   N       = size(xmean,1);
   I       = eye(N);
   popsize = floor(myeval(Options.popSize)*myeval(Options.incPopSize)^(irun-1)); 
   lambda  = popsize;

   % Adjust sigma0
   sigma0 = myeval(sigma0);
   if any(sigma0 <= 0)
      error('cmaes:sigma','Initial search volume (SIGMA) must be > 0 ');
   end
   if max(sigma0)/min(sigma0) > 1e+6
      error('cmaes:sigma','Initial search volume (SIGMA) badly conditioned');
   end
   if all(size(sigma0) == [N,2])
      sigma0 = 0.5*(sigma0(:,2) - sigma0(:,1));
   end

   % Evaluate options and set variable bounds
   Stop    = myeval(Options.Stop);
   maxdx   = myeval(Options.diffMaxChange);
   mindx   = myeval(Options.diffMinChange);
   lbounds = myeval(Options.lBounds);
   ubounds = myeval(Options.uBounds);
   if length(lbounds) == 1, lbounds = repmat(lbounds,N,1); end
   if length(ubounds) == 1, ubounds = repmat(ubounds,N,1); end
   if any(lbounds >= ubounds)
      error('cmaes:bounds','lower bound < upper bound must hold');
   end

   % Last possible setting of sigma0, in which case, tolerances are revaluated 
   if isempty(sigma0)
      if all(isfinite([lbounds,ubounds]))
         sigma0      = 0.3*(ubounds - lbounds);
         Stop.tolX   = myeval(Options.Stop.tolX);
         Stop.tolUpX = myeval(Options.Stop.tolUpX);
      else
         error('cmaes:sigma','Initial step sizes (SIGMA) not determined');
      end
   end

   % Error check all vector sizes
   if ~all(size(xmean) == [N,1])
      error('cmaes:sizes','Initial point has incompatible size');
   elseif ~(all(size(lbounds) == [N,1]) && all(size(ubounds) == [N,1]))
      error('cmaes:sizes','Variable bounds have incompatible size');
   elseif ~(all(size(sigma0) == [1,1])  || all(size(sigma0) == [N,1]))
      error('cmaes:sizes','Input parameter SIGMA has incompatible size');
   elseif size(Stop.tolX,2) > 1         || ~ismember(size(Stop.tolX,1), [1,N])
      error('cmaes:sizes','Option TolX has incompatible size');
   elseif size(Stop.tolUpX,2) > 1       || ~ismember(size(Stop.tolUpX,1), [1,N])
      error('cmaes:sizes','Option TolUpX has incompatible size');
   elseif size(maxdx,2) > 1             || ~ismember(size(maxdx,1), [1,N])
      error('cmaes:sizes','Option DiffMaxChange has incompatible size');
   elseif size(mindx,2) > 1             || ~ismember(size(mindx,1), [1,N])
      error('cmaes:sizes','Option DiffMinChange has incompatible size');
   end

   % Selection: Set number of parents/points and weights for recombination
   mu = myeval(Options.nParents);
   switch Options.recombWeights
   case 'equal'
      weights = ones(mu,1);
   case 'linear'
      weights = (mu:-1:1)';
   case 'superlinear'
      weights = log(mu+1) - log(1:mu)';
   otherwise
      error('cmaes:weights',['Recombination weights to be "' ...
            Options.recombWeights '" is not implemented']);
   end

   % Variance-effective size of mu
   mueff = sum(weights)^2/sum(weights.^2);
   if mueff == lambda
      error('cmaes:weights',['Combination of values for PopSize, ',...
            'ParentNumber, and RecombinationWeights is not reasonable']);
   end

   % Adaptation: Set time constants for covariance and step size control
   cc    = 4/(N+4);
   cs    = (mueff+2)/(N+mueff+3);
   mucov = mueff;
   ccov  = 2/(N+1.41)^2/mucov+(1-1/mucov)*min(1,(2*mueff-1)/((N+2)^2+mueff)); 

   % Set damping for controlling sigma (usually ~1, limits increase)
   damps = (1 + 2*max(0,sqrt((mueff-1)/(N+1))-1)) ...
         * max(0.3,1 - N/min(Stop.maxIter,Stop.maxFunEval/lambda)) + cs; 

   % Initialize dynamic internal state parameters
   pc    = zeros(N,1);
   ps    = zeros(N,1);
   sigma = max(sigma0);
   if length(sigma0) == 1
      sigma0 = sigma0*ones(N,1);
   end
   B  = I;
   D  = diag(sigma0/max(sigma0));
   BD = B*D;
   C  = BD*(BD)';
   fitness.hist    = nan(1,10+ceil(3*10*N/lambda));
   fitness.histsel = nan(1,10+ceil(3*10*N/lambda));

   % Initialize boundary handling
   bnd.isbounded  = (lbounds > -Inf) | (ubounds < Inf);
   bnd.exists     = any(bnd.isbounded);
   if bnd.exists
      xmean       = xintobounds(xmean,lbounds,ubounds);   % just in case
      bnd.weights = zeros(N,1);

      % Scaling is good, better in axis-parallel case, worse in rotated
      bnd.flgscale = 0;     % scaling will be omitted if zero 
      if bnd.flgscale
         bnd.scale = diag(C)/mean(diag(C));
      else
         bnd.scale = ones(N,1);
      end

      maxdx = min(maxdx,(ubounds - lbounds)/2);
      if any(sigma*sqrt(diag(C)) > maxdx)
         sigma = min(maxdx./sqrt(diag(C)));
      end
      dd  = diag(C);
      idx = (lbounds > -Inf) & (ubounds < Inf);
      if any(5*sigma*sqrt(dd(idx)) < ubounds(idx) - lbounds(idx))
         warning('cmaes:sigma', ['A SIGMA value is much smaller than its ', ...
                 'given bound interval.  For reasonable performance, ', ...
                 'values should lie between 0.2 and 0.5 of its bounds.']);
      end

      % delta fit for setting weights
      bnd.dfithist    = 1;
      bnd.validfitval = 0;
      bnd.iniphase    = 1;
   end

   % ooo initial feval, for output only
   if irun == 1 
      out.bestever.x     = xmean;
      out.bestever.f     = Inf;
      out.bestever.evals = nFunc;
   end
   if Options.evalInitialX
      fitness.hist(1)    = feval(f,xmean,varargin{:}); 
      fitness.histsel(1) = fitness.hist(1);
      nFunc              = nFunc + 1;
      if fitness.hist(1) < out.bestever.f 
         out.bestever.x     = xmean;
         out.bestever.f     = fitness.hist(1);
         out.bestever.evals = nFunc;
      end
      out.Cache.x = [out.Cache.x,xmean];
      out.Cache.f = [out.Cache.f,fitness.hist(1)];
   else
      fitness.hist(1)    = NaN; 
      fitness.histsel(1) = NaN; 
   end

   % Initialize random number generator
%   randn('state',myeval(Options.seed));
   rng(myeval(Options.seed),'v5normal');

   % Initialize further constants (expectation of ||N(0,I)||=norm(randn(N,1))
   chiN = N^0.5*(1 - 1/(4*N) + 1/(21*N^2));

   % Normalize recombination weights
   weights = weights/sum(weights);

   % -------------------- Generation Loop --------------------------------
   stopflag = zeros(length(stoplabels),1);
   nIter = 0;
   while ~any(stopflag)
      nIter = nIter + 1;
      if Options.debug >= 3, disp(['   CMA-ES ',int2str(nIter)]); end

      % Generate lambda offspring
      z = randn(N,lambda);
      x = repmat(xmean,1,lambda) + sigma*(BD*z);
      xvalid = xintobounds(x,lbounds,ubounds);

      % non-parallel evaluation of lambda offspring
      fitness.raw = nan(1,lambda);
      for k = 1:lambda
         while isnan(fitness.raw(k))
            fitness.raw(k) = feval(f,xvalid(:,k),varargin{:}); 
            out.Cache.x    = [out.Cache.x, xvalid(:,k)];
            out.Cache.f    = [out.Cache.f, fitness.raw(k)];
            if isnan(fitness.raw(k))
               z(:,k)      = randn(N,1);
               x(:,k)      = xmean + sigma*(BD*z(:,k));
               xvalid(:,k) = xintobounds(x(:,k),lbounds,ubounds);
            end
         end
      end
      nFunc = nFunc + lambda; % retries due to NaN are not counted
      fitness.sel = fitness.raw; 

      % Begin handling of variable bounds
      if bnd.exists

         % Get delta fitness values.  More precise: exp(mean(log(diag(C))))
         val = myprctile(fitness.raw, [25,75]);
         val = (val(2) - val(1)) / N / mean(diag(C)) / sigma^2;

         % Catch non-sensible values 
         if ~isfinite(val)
            warning('cmaes:fitness','Non-finite fitness range');
            val = max(bnd.dfithist);  
         elseif val == 0 % happens if all points are out of bounds
            val = min(bnd.dfithist(bnd.dfithist > 0)); 
         elseif bnd.validfitval == 0 % first sensible val
            bnd.dfithist    = [];
            bnd.validfitval = 1;
         end

         % Store delta fitness values
         bnd.maxhist  = length(bnd.dfithist) >= 20+(3*N)/lambda;
         bnd.dfithist = [bnd.dfithist(bnd.maxhist+1:end), val];
         [tx,ti]      = xintobounds(xmean,lbounds,ubounds);

         % Set initial weights
         if any(ti)
            if bnd.iniphase
               bnd.weights(bnd.isbounded) = 2.0002*median(bnd.dfithist);
               if ~bnd.flgscale
                  dd  = diag(C);
                  idx = find(bnd.isbounded);
                  dd  = dd(idx)/mean(dd); %  remove mean scaling
                  bnd.weights(idx) = bnd.weights(idx)./dd;
               end
               if bnd.validfitval && nIter > 2
                  bnd.iniphase = 0;
               end
            end

            % Increase weights for any variable that violates bounds 
            %  (and only if % xmean is moving away)
            tx  = xmean - tx; % distance from xmean to boundary
            idx = (ti~=0 & (sign(tx) == sign(xmean - xold)) ...
                  & abs(tx) > 3*max(1,sqrt(N)/mueff)*sigma*sqrt(diag(C)));
            if ~isempty(idx)
               bnd.weights(idx) = 1.2^(max(1, mueff/10/N))*bnd.weights(idx);
            end
         end

         % Calculate scaling biased to unity, product is one
         if bnd.flgscale
            bnd.scale = exp(0.9*(log(diag(C)) - mean(log(diag(C))))); 
         end

         % Assigned penalized fitness
         bnd.penalty = (bnd.weights./bnd.scale)'*(xvalid-x).^2; 
         fitness.sel = fitness.raw + bnd.penalty;
      end

      % Sort by fitness value and record short history of best fitness values
      [fitness.raw,fitness.idx]    = sort(fitness.raw);
      [fitness.sel,fitness.idxsel] = sort(fitness.sel);
      fitness.hist(2:end)          = fitness.hist(1:end-1);
      fitness.hist(1)              = fitness.raw(1);
      fitness.histsel(2:end)       = fitness.histsel(1:end-1);
      fitness.histsel(1)           = fitness.sel(1);

      % Calculate new xmean; this is selection and recombination 
      xold  = xmean;                           % for speed up of Eq. (2) and (3)
      xmean = x(:,fitness.idxsel(1:mu))*weights; 
      zmean = z(:,fitness.idxsel(1:mu))*weights; %==D^-1*B'*(xmean-xold)/sigma

      % Cumulation: update evolution paths
      ps   = (1-cs)*ps+(sqrt(cs*(2-cs)*mueff))*(B*zmean);               % Eq (4)
      hsig = norm(ps)/sqrt(1-(1-cs)^(2*nIter))/chiN < 1.4 + 2/(N+1);
%      pc   = (1-cc)*pc+hsig*(sqrt(cc*(2-cc)*mueff)/sigma)*(xmean-xold); % Eq (2)
      pc   = (1-cc)*pc+hsig*sqrt(cc*(2-cc)*mueff)*(B*D*zmean); % Eq (2)

      % Adapt covariance matrix: old matrix + rank-1 update + rank-mu update
      if ccov > 0                                                       % Eq (3)
         C = (1-ccov+(1-hsig)*ccov*cc*(2-cc)/mucov)*C ...
           + ccov*(1/mucov)*(pc*pc') + ccov*(1-1/mucov) ...
           * (x(:,fitness.idxsel(1:mu))-repmat(xold,1,mu))/sigma^2 ...
           * diag(weights)*(x(:,fitness.idxsel(1:mu))-repmat(xold,1,mu))';
      end

      % If ps is large and fitness is getting worse, remove momentum in ps.
      % This should rarely happen and is questionable in dynamic environments.
      if sum(ps.^2)/N > 1.5 + 10*sqrt(2/N) ...
                      && fitness.histsel(1) > max(fitness.histsel(2:3))
         ps = ps*sqrt(N*(1+max(0,log(sum(ps.^2)/N)))/sum(ps.^2));
      end

      % Adapt sigma
      sigma = sigma*exp((norm(ps)/chiN - 1)*cs/damps);                  % Eq (5)

      % Update B and D from C
      if ccov > 0 && mod(nIter,1/ccov/N/10) < 1

         % Get eigenvalues, normalized eigenvectors from symmetrized C
         C = triu(C) + triu(C,1)';
         [B,D] = eig(C);
         d = diag(D);
         if any(~isfinite(d)) || any(any(~isfinite(diag(B))))
            error('cmaes:eig',['function eig returned nonfinite eigenvalues,'...
                  ' cond(C) = ', num2str(cond(C)) ]);
         end

         % limit condition number of C to 1e+14 + 1
         if min(d) <= 0
            D(D < 0) = 0;
            temp = max(d)/1e+14;
            C = C + temp*I;
            D = D + temp*I; 
         end
         d = diag(D);
         if max(d) > 1e+14*min(d) 
            temp = max(d)/1e+14 - min(d);
            C = C + temp*I;
            D = D + temp*I;
         end

         % Store standard deviations in D
         D = diag(sqrt(diag(D)));
         % D = D / prod(diag(D))^(1/N);
         % C = C / prod(diag(D))^(2/N);
         BD = B*D;
      end

      % Adjust step size in the case of numerical precision problems
      c = diag(C);
      if any(sigma*sqrt(c) > maxdx)
         sigma = min(maxdx./sqrt(c));
      end
      while any(sigma*sqrt(c) < mindx)
         sigma = max(mindx./sqrt(c))*exp(0.05+cs/damps); 
      end
      temp = xmean == xmean + 0.2*sigma*sqrt(c);
      if any(temp)
         C = C + ccov*diag(c.*temp);
         sigma = sigma*exp(0.05+cs/damps); 
      end
      if all(xmean == xmean + 0.1*sigma*BD(:,1+floor(mod(nIter,N))))
         sigma = sigma*exp(0.2+cs/damps); 
      end
      if fitness.sel(1) == fitness.sel(1+ceil(0.1+lambda/4))
         sigma = sigma*exp(0.2+cs/damps); 
      end
      temp = [fitness.hist,fitness.sel(1)];
      if nIter > 2 && (max(temp) - min(temp)) == 0  
         sigma = sigma*exp(0.2+cs/damps); 
      end

      % Keep overall best solution
      if fitness.hist(1) < out.bestever.f
         out.bestever.x     = xvalid(:, fitness.idx(1));
         out.bestever.f     = fitness.hist(1);
         out.bestever.evals = nFunc + fitness.idx(1) - lambda;
      end

      % Set stop flags
      temp     = [fitness.sel,fitness.hist];
      stopflag = [sigma*max(diag(D)) == 0;   fitness.raw(1) <= Stop.fitness; 
                  nFunc >= Stop.maxFunEval; nIter >= Stop.maxIter;
                  all(sigma*(max(abs(pc), sqrt(diag(C)))) < Stop.tolX);
                  any(sigma*sqrt(diag(C)) > Stop.tolUpX);
                  nIter > 2 && (max(temp) - min(temp)) < Stop.tolFun;
                  exist('ABORT.txt','file')];
   end % while, end generation loop

   % Output best points
   out.stopflag(irun) = stoplabels(find(stopflag)); %#ok
   xmin = xvalid(:,fitness.idx(1));
   fmin = fitness.raw(1);

   % Stop run immediately for 'manual' or 'maxfunevals' condition is reached
   if any([stopflag(3),stopflag(7)])
      break
   end
end
return

%===============================================================================
% xintobounds:  Enforce lower and upper bounds on the variables stored in x.
% ------------------------------------------------------------------------------
% VARIABLES
%  x       = column vector or matrix of variable values
%  idx     = indices of bound violations
%  lbounds = vector of lower bounds
%  ubounds = vector of upper bounds
%  idx1    = indices of lower bound violations
%  idx2    = indices of upper bound violations
%===============================================================================
function [x,idx] = xintobounds(x,lbounds,ubounds)

lbounds = repmat(lbounds,1,size(x,2));
ubounds = repmat(ubounds,1,size(x,2));
idx1    = x < lbounds;
idx2    = x > ubounds;
x(idx1) = lbounds(idx1);
x(idx2) = ubounds(idx2);
idx     = idx2-idx1;
return

%===============================================================================
% myprctile: Computes the percentiles in vector perc from vector X returns
%            vector with length(res) == length(perc)
%===============================================================================
function res = myprctile(X, perc)

% Error check
if min(size(X)) > 1
   error('cmaes:myprctile','X must not be a matrix');
end
if min(size(perc)) > 1
   error('cmaes:myprctile','perc must not be a matrix');
end
if any(perc < 0) || any(perc > 100)
   error('cmaes:myprctile','perc must contain values in the range [0,100]');
end

% Process input variables
perc = perc(:);
sX   = sort(X);

% Convert percentiles
N         = length(X);
prctiles  = 100*((1:N)-0.5)/N;
res       = nan(N,length(perc));
res(perc <= prctiles(1))   = sX(1);
res(perc >= prctiles(end)) = sX(N);
ind       = find(perc > prctiles(1) & perc < prctiles(end));

% Find largest index smaller than requested percentile and linearly interpolate
for k = 1:length(ind)
   i = max(find(perc(ind(k)) > prctiles)); %#ok
   res(k) = sX(i) + (sX(i+1)-sX(i))*(perc(ind(k))-prctiles(i)) ...
                    /(prctiles(i+1)-prctiles(i));
end
return
