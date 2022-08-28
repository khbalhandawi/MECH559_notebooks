function Defaults = mads_defaults(varargin)
%MADS_DEFAULTS  Set default values for MADS parameters and options.
%
%   Syntax:
%      DEFAULTS = mads_defaults(TYPEPROBLEM)
%
%   Description:
%      MADS_DEFAULTS assigns default values for all variables passed into the
%      MADS optimizer.  It stores these values in a structure named DEFAULTS.
%      This function is called by several different NOMADm functions.
%      TYPEPROBLEM is a string that is set to either "Truth" or "Surrogate",
%      the former being used in almost all cases.  The only time "Surrogate" is
%      used is when the MADS optimizer is called within its own search function
%      to optimize a surrogate problem.
%
%   Note:
%      It is strongly recommended that these variables not be edited, as it is
%      the only place where these defaults are set.  This function is not meant
%      to be run independently.
%
%   See also MADS, MADS_BATCH, NOMADM

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
% mads_Defaults:  Sets default values for most of the MADS parameters
% ------------------------------------------------------------------------------
% Called by: NOMADm, nomadm_functions, mads, madsBatch
% VARIABLES:
%  Defaults            = structure containing all the default values
%  Labels              = long text labels for Search and Poll parameters
%    .file             =   labels for types of user files
%    .scale            =   labels for types of mesh direction scaling
%    .sensor           =   labels for sensor dimension in sensor placement
%    .filter           =   labels for types of filters
%    .degeneracyScheme =   labels for degenerate linear constraints schemes
%    .search           =   labels for types of Searches
%    .daceRegression   =   labels for types of DACE regression functions
%    .daceCorrelation  =   labels for types of DACE correlation functions
%    .nwKernel         =   labels for types of Nadaraya-Watson kernel functions
%    .rbfKernel        =   labels for types of radial basis function kernels
%    .rbfPoly          =   labels for types of radial basis function polynomials
%    .poll             =   labels for types of Poll strategies
%    .pollOrder        =   labels for types of Poll order strategies
%    .pollCenter       =   labels for types of Poll centers
%    .Parameters       =   labels for MADS parameters
%      .mesh           =     labels for mesh parameters
%      .term           =     labels for termination criteria parameters
%      .other          =     labels for other MADS parameters
%      .Stoch          =     labels for stochastic parameters
%  Types               = short text codes for Search and Poll choices
%    .file             =   types of user files
%    .degeneracyScheme =   types of schemes for degenerate linear constraints
%    .search           =   types of Searches
%    .daceRegression   =   types of DACE regression functions
%    .daceCorrelation  =   types of DACE correlation functions
%    .nwKernel         =   types of Nadaraya-Watson kernel functions
%    .pollStrategy     =   types of Poll strategies
%    .pollOrder        =   types of Poll order strategies
%  FileExt             = string file name suffixes
%    .F                =   functions file suffix
%    .X                =   closed constraints file
%    .O                =   Omega file suffix
%    .I                =   initial points file suffix
%    .N                =   discrete neighbor file suffix
%    .P                =   user parameter file suffix
%    .C                =   cache file suffix
%    .S                =   session file suffix
%    .D                =   diary or debug suffix
%  maxSearches         = maximum number of different Search types in a run
%  Choice              = integer choices of allowable strategies
%    .search           =   choice of Search type 
%    .daceRegr         =   choice of DACE regression function
%    .daceCorr         =   choice of DACE correlation function
%    .Poll.type        =   choice of Poll strategy
%    .Poll.order       =   choice of Poll order strategy
%    .Poll.center      =   choice of Poll center
%    .degeneracyScheme =   choice of scheme for degenerate linear constraints
%  Options             = values for the MADS Options structure
%    .nSearches        =   number of different Search types used
%    .Search           =   Search parameters
%    .dace             =   DACE parameters
%      .reg            =     regression function handle 
%      .corr           =     correlation function handle
%      .theta          =     initial guess of the theta fitting parameter
%      .lower          =     lower bound for the theta parameter
%      .upper          =     upper bound for the theta parameter
%      .isotropic      =     flag indicating isotropic theta values
%    .nw               =   NW parameters
%      .kernel         =     string name of NW kernel function 
%      .sigma          =     initial guess of the sigma fitting parameter
%      .lower          =     lower bound for the sigma parameter
%      .upper          =     upper bound for the sigma parameter
%    .rbf              =   RBF parameters
%      .kernel         =     string name of RBF kernel function
%      .poly           =     string name of RBF polynomial type
%    .Poll.type        =   short text code for default Poll strategy
%    .Poll.order       =   short text code for default Poll order
%    .Poll.center      =   numeric code for default choice of Poll center
%    .Poll.complete    =   turns on/off complete Polling
%    .EPoll.completeN  =   turns on/off complete discrete neighbor Polling
%    .EPoll.complete   =   turns on/off complete Extended Polling
%    .EPoll.fTrigger   =   f-value Extended Poll trigger
%    .EPoll.hTrigger   =   h-value Extended Poll trigger
%    .SurOptimizer     =   string identifying the surrogate optimizer
%    .Term             =   values for termination criteria
%      .delta          =     minimum mesh size
%      .nIter          =     maximum number of iterations
%      .nFunc          =     maximum number of function evaluations
%      .time           =     maximum CPU time
%      .nFails         =     maximum number of consecutive Poll failures
%    .TermFlag         =   flags to turn on/off termination criteria
%      .nIter          =     turns on/off number of iterations
%      .nFunc          =     turns on/off number of function evaluations
%      .time           =     turns on/off CPU time
%      .nFails         =     turns on/off number of consecutive Poll failures
%    .loadCache        =   turns on/off loading of a pre-existing Cache
%    .countCache       =   flag for counting Cache points as function calls
%    .countEval        =   flag for counting X-infeasible as function evals
%    .useFilter        =   turns on/off filter for nonlinear constraints
%    .degeneracyScheme =   scheme for handling degenerate linear constraints
%    .removeRedundancy =   removes redundant linear constraints
%    .runSensor        =   sensor dimension, if a sensor placement problem
%    .runStochastic    =   flag for running as stochastic optimization problem
%    .accelerate       =   flag for accelerating mesh refinement
%    .scale            =   base for logarithmic scaling (0 = no scaling)
%    .plotFilter       =   turns on/off real-time filter plot
%    .plotHistory      =   turns on/off history and real-time history plot
%    .delta0           =   initial mesh size
%    .deltaMax         =   maximum mesh size
%    .meshRefine       =   mesh refinement factor
%    .meshCoarsen      =   mesh coarsening factor
%    .tolCache         =   tolerance for identifying points already in Cache
%    .hmin             =   minimum h-value of an infeasible point
%    .hmax             =   maximum h-value of a filter point
%    .runOneIteration  =   flag for running MADS one iteration at a time
%    .runUntilFeasible =   flag for running MADS only unitl feasible
%    .runCount         =   run counter
%    .TermRel          =   flag for multiplying Term.delta by initial delta
%    .tolBind          =   tolerance for ID-ing active linear constraints
%    .nRuns            =   number of runs for stochastic problem
%    .hplothandle      =   handle for history plot axes
%    .fplothandle      =   handle for filter plot axes
%    .stophandle       =   handle for Stop Run pushbutton
%  Sur                 = structure with the same fields as Defaults.Options
%  typeProblem         = string ID: "Truth" or "Surrogate"
%*******************************************************************************

% Labels for Search, Poll, and MADS Parameters
Labels.version          =  '5.0';
Labels.file             = {'Functions'; 'Initial Points'; 'Closed Constraints';
                           'Linear Constraints'; 'Neighbors'; 'Parameter'};
Labels.scale            = {'No Scaling';
                           'Logarithmic Base-2 Scaling';
                           'Logarithmic Base-10 Scaling'};
Labels.sensor           = {'Off';'1-D';'2-D';'3-D'};
Labels.debug            = {'Off';'Level 1';'Level 2';'Level 3'};
Labels.plotHistory      = {'Off';'On';'Plot Real-Time'};
Labels.degeneracyScheme = {'Sequential Selection';
                           'Random Selection';
                           'Proximity to Constraint Boundary';
                           'Full Enumeration'};
Labels.filter           = {'No Filter';
                           'Traditional Filter';
                           'Progressive Barrier'};
Labels.multiFi          = {'No Multifidelity Algorithm';
                           'Apply Recursive MADS Algorithm';
                           'Add Fidelity Level as Categorical Variable '};
Labels.search           = {'None';
                           'Central Composite Design Sampling';
                           'Latin Hypercube Sampling';
                           'Grid Sampling on coarse mesh';
                           'Particle Swarm';
                           'CMAES';
                           'Standard Poll about N filter points';
                           'Gradient Poll about N filter points';
                           'Kriging Surrogate';
                           'Radial Basis Function Surrogate';
                           'Nadaraya-Watson Surrogate';
                           'SPS with Kriging Surrogate';
                           'SPS with Radial Basis Function Surrogate';
                           'SPS with Nadaraya-Watson Surrogate';
                           'Custom Search Function';
                           'Custom Surrogate Function'};
Labels.cmaesRecombWeights = {'superlinear';'linear';'equal'};
Labels.daceRegression   = {'Zero-Order Polynomial';
                           'First-Order Polynomial';
                           'Second-Order Polynomial'};
Labels.daceCorrelation  = {'Exponential';
                           'General Exponential';
                           'Gaussian';
                           'Local Support, Linear';
                           'Local Support, Spherical';
                           'Local Support, Cubic Polynomial';
                           'Local Support, Cubic Spline'};
Labels.mleOpt           = {'Matlab MultiStart';
                           'Matlab GlobalSearch';
                           'Matlab fmincon';
                           'MADS';
                           'DACE boxmin';
                           'Grid Search';
                           'Grid Search with MultiStart'};
Labels.nwKernel         = {'Gaussian';
                           'Uniform';
                           'Triangle';
                           'Epanechnikov';
                           'Quartic';
                           'Tri-weight';
                           'Cosinus'};
Labels.rbfKernel        = {'Bi-harmonic';
                           'Tri-harmonic';
                           'Gaussian';
                           'Multiquadric';
                           'Inverse Multiquadric';
                           'Thin Plate Spline'};
Labels.rbfPoly          = {'Zero-Order Polynomial';
                           'Second-Order Polynomial';
                           'Reduced Second-Order Polynomial';
                           'Third-Order Polynomial';
                           'Reduced Third-Order Polynomial'};
Labels.poll             = {'Standard 2n e-directions';
                           'Standard n+1 directions';
                           'Custom 2n directions';
                           'Custom n+1 directions';
                           'MADS Random 2n directions';
                           'MADS Random n+1 directions';
                           'MADS Random 2n orthogonal directions';
                           'MADS Random n+1 orthogonal directions';
                           'MADS Random 2-point poll';
                           'MADS Deterministic 2n orthogonal directions';
                           'MADS Deterministic n+1 orthogonal directions';
                           'MADS Deterministic 2-point poll';
                           'Gradient-pruned 2n e-directions';
                           'Gradient-pruned n+1 directions';
                           'Gradient-pruned 3^n using L-1 gradient';
                           'Gradient-pruned 3^n using L-2 gradient';
                           'Gradient-pruned 3^n using L-Inf gradient';
                           'Gradient-pruned 3^n Descent Vectors'};
Labels.pollOrder        = {'Consecutive';
                           'Alternating';
                           'Random';
                           'Dynamic';
                           'Dynamic Ranked';
                           'Simplex Gradient';
                           'Surrogate Ranked';
                           'Custom'};
Labels.pollCenter       = {'Best Feasible Point';
                           'Least Infeasible Point';
                           '2nd Least Infeasible Filter Point';
                           '3rd Least Infeasible Filter Point';
                           '4th Least Infeasible Filter Point';
                           '5th Least Infeasible Filter Point'};
Labels.Parameters.term  = {'Convergence Tolerance (Mesh Size):';
                           'Maximum Number of Iterations:';
                           'Maximum Number of Function Calls:';
                           'Maximum CPU Time:';
                           'Maximum Number of Consecutive Poll Failures:'};
Labels.Parameters.mesh  = {'Initial Mesh Size:';
                           'Maximum Mesh Size:';
                           'Mesh Refinement Factor:';
                           'Mesh Coarsening Factor:';
                           'Cache Tolerance:'};
Labels.Parameters.other = {'Minimum Filter Constraint Violation:';
                           'Maximum Filter Constraint Violation:';
                           'MADS-PB Frame Center Trigger:';
                           'MVP Objective Extended Poll Trigger:';
                           'MVP Constraints Extended Poll Trigger:';
                           'Number of Runs (for stochastic problems):'};
Labels.coolCats         = {'Charles Audet', 'Keith Berrier','Olga Brezhneva',...
                           'Gilles Couture','Ana Custodio', 'John Dennis',   ...
                           'Thierry Dalon', 'John Dunlap', ...
                           'Aran Garcia-Lekue','Alison Marsden', ...
                           'Todd Paciencia','Rachael Pingel Robison', ...
                           'Jacob Sondergaard','Todd Sriver','Luis Vicente'};

% Surrogate optimizer labels
Labels.optimizer        = {'FMINCON (MATLAB Optimization Toolbox)';
                           'MADS';
                           'Particle Swarm Optimizer';
                           'CMA-ES (Evolutionary Strategy)';
                           'Custom'};
if ~license('test','Optimization_Toolbox'), Labels.optimizer(1) = []; end

% Names of help files 
HelpDoc = struct('nomadm', 'nomadm_help.pdf', ...
                 'dace',   'dace.pdf',        ...
                 'nw',     'nw_help.pdf',     ...
                 'rbf',    'rbf_help.pdf',    ...
                 'changes','nomadm.txt',      ...
                 'license','gpl.txt');

% Short text codes for Search, DACE Surrogates, and Poll parameters
Types.file      = {'','_x0','_X','_Omega','_N','_Param'};
Types.search    = {'None','CCD','LHS','Mesh','PS','CMAES','SPollI','GPollI', ...
                   'DACE','RBF','NW','SPS-DACE','SPS-RBF','SPS-NW',...
                   'Custom','CustomS'};
Types.daceRegr  = {'regpoly0','regpoly1','regpoly2'};
Types.daceCorr  = {'correxp','correxpg','corrgauss','corrlin',...
                   'corrspherical','corrcubic','corrspline'};
Types.mleOpt    = {'MultiStart','GlobalSearch','fmincon','mads','boxmin', ...
                   'grid','gridOpt'};
Types.rbfPoly   = {'regpoly0','regpoly2','regpoly2reduced', ...
                   'regpoly3','regpoly3reduced'};
Types.rbfKernel = {'bi-harmonic','tri-harmonic','gaussian','multiquadric',...
                   'inv-multiquadric','thinplatespline'};
Types.nwKernel  = {'gaussian','uniform','triangle','Epanechnikov','quartic', ...
                   'triweight','cosinus'};
Types.cmaesRecombWeights = {'superlinear','linear','equal'};
Types.cmaesStopLabels = {};
Types.optimizer = {'fmincon','mads','ps','cmaes','custom'};
Types.poll      = {'Standard_2n',     'Standard_n+1',              ...
                   'Custom_2n',       'Custom_n+1',                ...
                   'MADS_2n',         'MADS_n+1',                  ...
                   'OrthoMADSr_2n',   'OrthoMADSr_n+1', 'MADSr_2', ...
                   'OrthoMADS_2n',    'OrthoMADS_n+1',  'MADS_2',  ...
                   'Gradient_2n',     'Gradient_n+1',              ...
                   'Gradient_3n_L1',  'Gradient_3n_L2',            ...
                   'Gradient_3n_LInf','Gradient_3n2n'};
Types.pollOrder = {'Consecutive',  'Alternating',    'Random',   'Dynamic', ...
                   'DynamicRanked','SimplexGradient','Surrogate','Custom'};
Types.degeneracyScheme = {'sequential','random','closest','full'};
Types.plotColors       = 'kbrmgkbrmg';

% File Extentions
FileExt.F = '';
FileExt.I = '_x0';
FileExt.X = '_X';
FileExt.O = '_Omega';
FileExt.N = '_N';
FileExt.P = '_Param';
FileExt.C = '_Cache.mat';
FileExt.S = '_Session.mat';
FileExt.H = '_History';
FileExt.D = '_Debug.txt';

% Choices for Search and Poll
maxSearches           = 8;
Options.nSearches     = 2;
Choice.optimizer      = 1;
Choice.Poll.type      = 1;
Choice.Poll.order     = 1;
Choice.Poll.center    = 1;
Choice.EPoll          = Choice.Poll;
Options.SurOptimizer  = Types.optimizer{Choice.optimizer};
Options.mvp1Surrogate = 1;
for k = 1:maxSearches
   Choice.search(k)                  = 1;
   Choice.daceRegr(k)                = 1;
   Choice.daceCorr(k)                = 3;
   Choice.mleOpt(k)                  = 1;
   Choice.nwKernel(k)                = 1;
   Choice.rbfKernel(k)               = 1;
   Choice.rbfPoly(k)                 = 1;
   Choice.cmaesRecombWeights(k)      = 1;
   Options.Search(k).type            = Types.search{Choice.search(k)};
   Options.Search(k).label           = Labels.search{Choice.search(k)};
   Options.Search(k).nIter           = 1;
   Options.Search(k).nPoints         = 1;
   Options.Search(k).complete        = 1;
   Options.Search(k).sfile           = {'',''};
   Options.Search(k).file            = '';
   Options.Search(k).local           = 0;
   Options.Search(k).merit           = 0;
   Options.Search(k).param           = 0;
   Options.Search(k).cblgs.nGoal     = 5;
   Options.Search(k).cblgs.nCloud    = 5000;
   Options.Search(k).cblgs.maxCloud  = 100000;
   Options.Search(k).cblgs.tolFeas   = 0.01;
   Options.Search(k).cblgs.maxIter   = 20;
   Options.dace(k).reg               = 'regpoly0';
   Options.dace(k).corr              = 'corrgauss';
   Options.dace(k).maxCondR          = exp(16);
   Options.dace(k).minNormDL         = 1e-3;
   Options.dace(k).mleTolBounds      = 1e-8;
   Options.dace(k).mleTolFinDiff     = 1e-1;
   Options.dace(k).mleOpt.type       = Types.mleOpt{Choice.mleOpt(k)};
   Options.dace(k).mleOpt.nPoints    = 10;
   Options.dace(k).isotropic         = 1;
   Options.nw(k).kernel              = Types.nwKernel{Choice.nwKernel(k)};
   Options.nw(k).sigma               = 0.5;
   Options.nw(k).lower               = 0.1;
   Options.nw(k).upper               = 3;
   Options.rbf(k).kernel             = Types.rbfKernel{Choice.rbfKernel(k)};
   Options.rbf(k).poly               = Types.rbfPoly{Choice.rbfPoly(k)};
   Options.cmaes(k).sigma            = 1;
   Options.cmaes(k).Stop.fitness     = -Inf;
   Options.cmaes(k).Stop.maxFunEval  =  Inf;
   Options.cmaes(k).Stop.maxIter     = '1e+3*(N+5)^2/sqrt(popsize)';
   Options.cmaes(k).Stop.tolX        = '1e-12*max(sigma0)';
   Options.cmaes(k).Stop.tolUpX      = '1e8*max(sigma0)';
   Options.cmaes(k).Stop.tolFun      = 1e-12;
   Options.cmaes(k).lBounds          = -Inf;
   Options.cmaes(k).uBounds          =  Inf;
   Options.cmaes(k).diffMaxChange    =  Inf;
   Options.cmaes(k).diffMinChange    = 0;
   Options.cmaes(k).evalInitialX     = 1;
   Options.cmaes(k).incPopSize       = 2;
   Options.cmaes(k).popSize          = '(4 + floor(3*log(N)))';
   Options.cmaes(k).nParents         = 'floor(popsize/2)';
   Options.cmaes(k).nRestarts        = 0;
   Options.cmaes(k).recombWeights    = 'superlinear';
   Options.cmaes(k).science          = 0;
   Options.cmaes(k).seed             = 'sum(100*clock)';
end
Options.Search(Options.nSearches).nIter = Inf;
Choice.degeneracyScheme    = 1;
Options.Poll.type          = Types.poll{Choice.Poll.type};
Options.Poll.label         = Labels.poll{Choice.Poll.type};
Options.Poll.complete      = 0;
Options.Poll.order         = Types.pollOrder{Choice.Poll.order};
Options.Poll.center        = Choice.Poll.center - 1;
Options.EPoll              = Options.Poll;
Options.EPoll.completeN    = 0;
Options.EPoll.fTrigger     = 0.01;
Options.EPoll.hTrigger     = 0.05;
Options.degeneracyScheme   = Types.degeneracyScheme{Choice.degeneracyScheme};
Options.lhsStrength        = 2;
Options.lhsBinFactor       = 1;

% MADS Parameter Values
Options.tolList           = reshape([1; 0.5; 0.2]*(10.^(0:-1:-14)),3*15,1);
Options.Term.delta        = 1e-4;
Options.Term.nIter        = 1000;
Options.Term.nFunc        = 50000;
Options.Term.time         = 3600;
Options.Term.nFails       = 50;
Options.TermFlag.delta    = 1;
Options.TermFlag.nIter    = 0;
Options.TermFlag.nFunc    = 1;
Options.TermFlag.time     = 0;
Options.TermFlag.nFails   = 0;
Options.TermFlag.relative = 0;
Options.loadCache         = 1;
Options.countCache        = 1;
Options.countNotInX       = 0;   % Count X-infeasible points as function evals
Options.useFilter         = 1;   % Traditional multipoint filter
Options.useMultiFi        = 0;   % Multifidelity options
Options.removeRedundancy  = 0;
Options.dimSensor         = 0;   % Dimension of each sensor
Options.runStochastic     = 0;
Options.fixCategorical    = 0;
Options.accelerate        = 0;
Options.scale             = 2;   % Base of logarithmic scaling (0 = none)
Options.scaleVariables    = 0;
Options.saveHistory       = 0;
Options.plotHistory       = 1;
Options.plotFilter        = 1;
Options.plotColor         = Types.plotColors(1);
Options.delta0            = 1.0;
Options.deltaMax          = Inf;
Options.meshRefine        = 0.5;
Options.meshCoarsen       = 1.0;
Options.tolCache          = Options.Term.delta;
Options.hmin              = 0.0;
Options.hmax              = 1.0;
Options.hRho              = 0.01;
Options.nRuns             = 1;
Options.runOneIteration   = 0;
Options.runUntilFeasible  = 0;
Options.tolBind           = 0.05;
Options.debug             = 0;
Options.digits            = 10;
Options.seed              = 0;
Options.abortFile         = 'ABORT.txt';

% Plot handles
Options.hplothandle = [];
Options.fplothandle = [];
Options.stophandle  = [];

% Surrogate defaults
Sur.Options = Options;
Sur.Options.nSearches         = 0;
Sur.Options.Search            = [];
Sur.Options.Poll.type         = 'Standard_2n';
Sur.Options.Poll.label        = 'Standard 2n e-directions';
Sur.Options.Poll.order        = 'Consecutive';
Sur.Options.Poll.center       = 0;
Sur.Options.Poll.complete     = 0;
Sur.Options.EPoll             = Sur.Options.Poll;
Sur.Options.EPoll.completeN   = 0;
Sur.Options.EPoll.fTrigger    = 0.01;
Sur.Options.EPoll.hTrigger    = 0.05;
Sur.Options.Term.delta        = 1e-3;
Sur.Options.Term.nIter        = 50;
Sur.Options.Term.nFunc        = 5000;
Sur.Options.Term.time         = Inf;
Sur.Options.Term.nFails       = Inf;
Sur.Options.TermFlag.delta    = 1;
Sur.Options.TermFlag.nIter    = 1;
Sur.Options.TermFlag.nFunc    = 1;
Sur.Options.TermFlag.time     = 0;
Sur.Options.TermFlag.nFails   = 0;
Sur.Options.TermFlag.relative = 0;
Sur.Options.loadCache         = 0;
Sur.Options.countCache        = 0;
Sur.Options.plotFilter        = 0;
Sur.Options.plotHistory       = 0;
Sur.Options.delta0            = Options.delta0;
Sur.Options.deltaMax          = Options.deltaMax;
Sur.Options.meshRefine        = Options.meshRefine;
Sur.Options.meshCoarsen       = Options.meshCoarsen;
Sur.Options.nRuns             = 1;
Sur.Options.runOneIteration   = 0;
Sur.Options.runUntilFeasible  = 0;
Sur.Options.useFilter         = 0;
Sur.Options.tolBind           = sqrt(Sur.Options.Term.delta);
Sur.Options.tolCache          = Sur.Options.Term.delta;
Sur.Options.debug             = Options.debug;

% Construct Default structure
Defaults.Labels      = Labels;
Defaults.Types       = Types;
Defaults.FileExt     = FileExt;
Defaults.Choice      = Choice;
Defaults.HelpDoc     = HelpDoc;
Defaults.maxSearches = maxSearches;

% Determine if a "Truth" or "Surrogate" function is to be optimized 
if nargin > 0
   Defaults.typeProblem = varargin{1};
else
   Defaults.typeProblem = 'TRUTH';
end

% Assign either Truth or Surrogate default values
switch upper(Defaults.typeProblem)
case 'TRUTH'
   Defaults.Options   = Options;
   Defaults.nameCache = 'CACHE';
case 'SURROGATE'
   Defaults.Options   = Sur.Options;
   Defaults.nameCache = 'sCACHE';
otherwise
   error('mads:defaults', ...
         'Bad argument for mads_defaults: Must be "TRUTH" or "SURROGATE"');
end
Defaults.Options.Sur = Sur.Options;
return
