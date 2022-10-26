function nomadm
%NOMADM   Execute the NOMADm graphic user interface (GUI).
%
%   Syntax:
%      nomadm
%
%   Description:
%      NOMADM launches the NOMADm GUI, which controls the setup of an
%      optimization problem, setting of various algorithm parameters and
%      options, running of the MADS algorithm for numerically solving the
%      optimization problem, and viewing results.
%
%   See also MADS, MADS_BATCH, MADS_DEFAULTS, NOMADM_COMPILE

%*******************************************************************************
%   Copyright (c) 2001-2017 by Mark A. Abramson
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
%   Last modified, 17 June 2017
%
%   Author information:
%   Mark A. Abramson, LtCol (ret.), USAF, PhD
%   Utah Valley University, Orem, Utah USA
%   Abramson.Mark@gmail.com
%*******************************************************************************

%*******************************************************************************
% nomadm: Runs the NOMADm user interface and associated callback functions.
% ------------------------------------------------------------------------------
% CALLED FUNCTIONS (CB = callback function attached to GUI object):
%  nomadm_gui               = Launches the NOMADm GUI.
%    selectProblem          =   Choose an optimization problem (CB)
%    editFile               =   Edit a user file (CB)
%    quitGUI                =   Close the GUI and quit the program (CB)
%    selectPoll             =   Select Poll strategy, center, or order (CB)
%    loadSession            =   Load user options from file (CB)
%    resetSession           =   Reset parameters to default values (CB)
%    saveSession            =   Save user options to file (CB)
%    saveCache              =   Save the Cache to file for later use (CB)
%    deleteFile             =   Delete a Cache or Session file (CB)
%    clearRuns              =   Clears all data from with previous runs (CB)
%    getHelp                =   View a help file (CB)
%    getAbout               =   View information about this software (CB)
%    copyPlot               =   Copy a plot to its own window (BD)
%    runMADS                =   Run the MADS optimizer (CB)
%      loadMADS             =     Transfer GUI values to MADS variables
%    edit_gui               =   Launch the Edit screen (CB)
%      loadEdit             =     Updates screen for edited parameter choices
%    search_gui             =   Launch the Search screen GUI
%      loadUserSearchFile   =     Load a user-specified Search file (CB)
%      loadSearchOptions    =     Load user Search step choices (CB)
%        modifySearchScreen =       Change appearance of the Search screen (CB)
%    cmaes_gui              =   Launch the CMA-ES figure window (CB)
%      loadCMAESOptions     =     Load user CMA-ES options (CB)
%    dace_gui               =   Launch the DACE figure window
%    nw_gui                 =   Launch the NW figure window (CB)
%    rbf_gui                =   Launch the RBF figure window (CB)
%      loadSurrogateOptions =     Load user surrogate options (CB)
%    updateScreen           =   Update contents of the NOMADm figure window
%      updateSearchLabels   =     Update Search labels on the NOMADm window
%    results_gui            =   Launch Results GUI, view results of run (CB)
%      solution_gui         =     Launch Solution GUI, view solution (CB)
%        changeView         =       Change view page using Previous/Next buttons
% ------------------------------------------------------------------------------
% VARIABLES (only for the nomadm function):
%  h                  = handle of pre-existing NOMADm figure window
%  gui_var            = structure of all GUI variables
%    .path            =   Current path
%    .runCount        =   current MADS run number
%    .runMax          =   maximum allowed MADS run number
%    .noGrad          =   flag indicating initial availability of gradients
%    .Defaults        =   substructure of parameter default values
%    .Labels          =   substructure of NOMADm menu labels
%    .Types           =   substructure of NOMADm types 
%    .FileExt         =   substructure of possible file extensions
%    .HelpDoc         =   substructure of help file names
%    .maxSearches     =   maximum allowable number of Search types used
%    .nameCache       =   name of the Cache that is saved as appdata
%  guiFunction        = cell array of NOMADm function handle names
%  gui                = structure of GUI object handles
%    .fig             =   handle for the main figure window
%    .func            =   structure of handles for all NOMADm functions
%*******************************************************************************

% Before launching the NOMADm GUI, make sure it is not already running
warning('off','MATLAB:dispatcher:InexactMatch');
h = findobj(0,'type','Figure','Tag','NOMADm_GUI');
if ~isempty(h)
   disp('Only one NOMADm window may run at a time');
   figure(h);
   return
end
clc

% Declare basic GUI variables
gui_var.path        = pwd;
gui_var.runCount    = 0;
gui_var.runMax      = 10;
gui_var.noGrad      = 1;
if ~isdeployed, addpath(pwd); end

% Set up default values
gui_var.Defaults    = mads_defaults('TRUTH');
gui_var.Labels      = gui_var.Defaults.Labels;
gui_var.Types       = gui_var.Defaults.Types;
gui_var.FileExt     = gui_var.Defaults.FileExt;
gui_var.HelpDoc     = gui_var.Defaults.HelpDoc;
gui_var.maxSearches = gui_var.Defaults.maxSearches;
gui_var.nameCache   = gui_var.Defaults.nameCache;
gui_var.typeProblem = gui_var.Defaults.typeProblem;

% Set function handles for all nomadm functions used in the GUI
guiFunction = { ...
   'selectProblem','editFile','quitGUI','selectPoll', ...
   'loadSession','resetSession','saveSession','saveCache','deleteFile', ...
   'clearRuns','getHelp','getAbout','copyPlot','runMADS','loadMADS', ...
   'updateScreen','edit_gui','search_gui','poll_gui','cmaes_gui', ...
   'dace_gui','nw_gui','rbf_gui','results_gui','solution_gui', ...
   'loadEdit','loadSearchOptions','loadPollOptions', ...
   'loadCMAESOptions','loadSurrogateOptions','loadUserSearchFile',...
   'changeView','modifySearchScreen','updateSearchLabels','uitoggle'};
for k = 1:length(guiFunction)
   gui.func.(guiFunction{k}) = str2func(guiFunction{k});
end

% Launch the NOMADm GUI and store GUI handles and variables as appdata
gui = nomadm_gui(gui,gui_var);
setappdata(gui.fig,'gui',gui);
setappdata(gui.fig,'gui_var',gui_var);
resetSession(1,[]);
end

%*******************************************************************************
% nomadm_gui:  Launch the NOMADm GUI.
% ------------------------------------------------------------------------------
% Calls: mads_defaults
% VARIABLES:
%  gui_var           = structure of GUI variables (descriptions above)
%  gui               = structure of all GUI object handles
%    .func           =   structure of function handles with fields as names
%    .fig            =   handle for the main figure window
%    .TermFlag       =   structure of handles for termination checkboxes
%      .nIter        =     handle for number of iterations
%      .nFunc        =     handle for number of function evaluations
%      .time         =     handle for CPU time
%      .nFails       =     handle for number of consecutive Poll failures
%    .problem        =   handle for text listing optimization problem name
%    .searchLabel(k) =   handles for Searches
%    .pollStrategy   =   handle for Poll strategy
%    .pollCenter     =   handle for Poll center
%    .pollOrder      =   handle for Poll order type
%    .Term           =   structure of handles for termination criteria
%      .delta        =     handle for poll size tolerance
%      .nIter        =     handle for maximum number of MADS iterations
%      .nFunc        =     handle for maximum number of function evaluations
%      .time         =     handle for maximum CPU time
%      .nFails       =     handle for max number of consec Poll failures
%    .delta0         =   handle for initial poll size
%    .deltaMax       =   handle for maximum poll size
%    .meshRefine     =   handle for mesh refinement factor
%    .meshCoarsen    =   handle for mesh coarsening factor
%    .tolCache       =   handle for Cache tolerance
%    .hmin           =   handle for minimum infeasible h-value
%    .hmax           =   handle for maximum filter h-value
%    .ePollXiF       =   handle for f-value Extended Poll trigger
%    .ePollXiH       =   handle for h-value Extended Poll trigger
%    .nRuns          =   handle for number of runs for stochastic problem
%    .runStatus      =   handle for GUI figure window status bar
%    .axesHistory    =   handle for the History plot axes
%    .axesFilter     =   handle for the Filter plot axes
%    .stopRun        =   handle for the Stop Run pushbutton
%    .resumeRun      =   handle for the Resume Run pushbutton
%    .Menu           =   structure of handles for the main menu
%      .Problem      =     handle for the Problem menu
%      .MADS         =     handle for the MADS menu
%      .Options      =     handle for the Options menu
%      .Session      =     handle for the Session menu
%      .Cache        =     handle for the Cache menu
%      .Run          =     handle for the Run menu
%      .Results      =     handle for the Results menu
%      .Help         =     handle for the Help menu
%    .ProblemMenu    =   structure of handles for Problem menu items
%    .MADSMenu       =   structure of handles for MADS menu items
%    .OptionsMenu    =   structure of handles for Options menu items
%    .SessionMenu    =   structure of handles for Session menu items
%    .CacheMenu      =   structure of handles for Cache menu items
%    .RunMenu        =   structure of handles for Run menu items
%    .ResultsMenu    =   structure of handles for Results menu items
%    .HelpMenu       =   structure of handles for Help menu items
%   onoff            = cell array containing "on" and "off" strings
%
%   MANY other variables, primarily ones that are GUI.fig objects
%*******************************************************************************
function gui = nomadm_gui(gui,gui_var)
onoff = {'on','off'};

% Set up main figure window
gui.fig = figure(...
   'Name',                               'NOMADm Optimization Software',...
   'Tag',                                'NOMADm_GUI', ...
   'DefaultUIControlUnits',              'normalized', ...
   'DefaultUIControlFontUnits',          'normalized', ...
   'DefaultUIControlFontName',           'Helvetica',  ...
   'DefaultUIControlFontSize',            0.72, ...,
   'DefaultUIControlStyle',              'text', ...
   'DefaultUIControlHorizontalAlignment','left', ...
   'Units',                              'normalized', ...
   'Position',                           [0.05 0.1 .9 .8], ...
   'MenuBar',                            'none', ...
   'NumberTitle',                        'off', ...
   'Color',                              [.8, .8, .8], ...
   'CloseRequestFcn',                    {gui.func.quitGUI, gui_var});

% object position parameters
headerHeight    = .035;
textHeight      = .03;
panelLeft       = .02;
panelWidth      = .435;
labelLeft       = .025;
labelShortWidth = .13;
labelLongWidth  = .34;
fieldShortLeft  = labelShortWidth + .03;
fieldLongLeft   = labelLongWidth  + .03;
fieldShortWidth = .08;
fieldLongWidth  = .29;

% Set up GUI status bar and panels
gui.statusbar = uipanel('Parent',gui.fig,'Position',[0, 0, 1, .04], ...
                                         'BorderType','beveledout');
gui.runStatus = uicontrol(gui.statusbar, 'String','No problem selected', ...
                                         'Position', [.01 0 .98 .88]);
gui.ProblemPanel = uipanel('Parent',gui.fig, ...
                           'Position',[panelLeft, .945, panelWidth, .045], ...
                           'BorderType','beveledout');
gui.MADSPanel    = uipanel('Parent',gui.fig, ...
                           'Position',[panelLeft, .675, panelWidth, .27], ...
                           'BorderType','beveledout');
gui.TermPanel    = uipanel('Parent',gui.fig, ...
                           'Position',[panelLeft, .515, panelWidth, .16], ...
                           'BorderType','beveledout');
gui.MeshPanel    = uipanel('Parent',gui.fig, ...
                           'Position',[panelLeft, .355, panelWidth, .16], ...
                           'BorderType','beveledout');
gui.FilterPanel  = uipanel('Parent',gui.fig, ...
                           'Position',[panelLeft, .165, panelWidth, .19], ...
                           'BorderType','beveledout');

% Display of the NOMADm object labels
uicontrol(gui.fig, 'String','Optimization Problem: ', ...
      'FontWeight','bold', ...
      'Position',     [labelLeft, .95, .16, headerHeight], ...
      'ToolTipString','Name of the currently selected optimization problem');
uicontrol(gui.fig, 'String','MADS Parameter Settings', ...
      'FontWeight','bold','HorizontalAlignment','center', ...
      'Position',     [labelLeft, .905, panelWidth-.025, textHeight], ...
      'ToolTipString','Any of these values may be changed via the MADS Menu');

uicontrol(gui.fig, 'String','Initial Search:', ...
      'Position',     [labelLeft, .87, labelShortWidth, textHeight], ...
      'ToolTipString','The Strategy used in the initial MADS Search step');
uicontrol(gui.fig, 'String','Search:', ...
      'Position',     [labelLeft, .84, labelShortWidth, textHeight], ...
      'ToolTipString','The Strategy used in the MADS Search step');

uicontrol(gui.fig, 'String','Poll Directions:', ...
      'Position',     [labelLeft, .74, labelShortWidth, textHeight], ...
      'ToolTipString', ...
      'MADS Poll directions must positively span the problem domain');
uicontrol(gui.fig, 'String','Poll Order:', ...
      'Position',     [labelLeft, .71, labelShortWidth, textHeight], ...
      'ToolTipString','The order in which the MADS Poll set is evaluated');
uicontrol(gui.fig, 'String','Poll Center:', ...
      'Position',     [labelLeft, .68, labelShortWidth, textHeight], ...
      'ToolTipString','The point around which the MADS Poll step is performed');
gui.TermFlag.delta  = uicontrol(gui.fig, ...
      'String',       gui_var.Labels.Parameters.term{1}, ...
      'Style',        'checkbox', ...
      'Value',        1, ...
      'Position',     [labelLeft, .64, labelLongWidth, textHeight]);
gui.TermFlag.nIter  = uicontrol(gui.fig, ...
      'String',       gui_var.Labels.Parameters.term{2}, ...
      'Style',        'checkbox', ...
      'Position',     [labelLeft, .61, labelLongWidth, textHeight]);
gui.TermFlag.nFunc  = uicontrol(gui.fig, ...
      'String',       gui_var.Labels.Parameters.term{3}, ...
      'Style',        'checkbox', ...
      'Position',     [labelLeft, .58, labelLongWidth, textHeight]);
gui.TermFlag.time   = uicontrol(gui.fig, ...
      'String',       gui_var.Labels.Parameters.term{4}, ...
      'Style',        'checkbox', ...
      'Position',     [labelLeft, .55, labelLongWidth, textHeight]);
gui.TermFlag.nFails = uicontrol(gui.fig, ...
      'String',       gui_var.Labels.Parameters.term{5}, ...
      'Style',        'checkbox', ...
      'Position',     [labelLeft, .52, labelLongWidth, textHeight]);
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.mesh{1}, ...
      'Position',     [labelLeft, .48, labelLongWidth, textHeight], ...
      'ToolTipString','The initial mesh size parameter value');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.mesh{2}, ...
      'Position',     [labelLeft, .45, labelLongWidth, textHeight], ...
      'ToolTipString','The maximum allowed mesh size');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.mesh{3}, ...
      'Position',     [labelLeft, .42, labelLongWidth, textHeight], ...
      'ToolTipString', ...
      'Mesh size is reduced by this factor when an iteration fails');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.mesh{4}, ...
      'Position',     [labelLeft, .39, labelLongWidth, textHeight], ...
      'ToolTipString', ...
      'Mesh size is increased by this factor when an iteration succeeds');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.mesh{5}, ...
      'Position',     [labelLeft, .36, labelLongWidth, textHeight], ...
      'ToolTipString', ...
      ['Any two points whose distance is less than this value are', ...
       ' assumed identical']);
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.other{1}, ...
      'Position',     [labelLeft, .32, labelLongWidth, textHeight], ...
      'ToolTipString', ...
      'Minimum constraint violation function value of an infeasible point');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.other{2}, ...
      'Position',     [labelLeft, .29, labelLongWidth, textHeight], ...
      'ToolTipString', ...
      'Maximum constraint violation function value of any filter point');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.other{3}, ...
      'Position',     [labelLeft, .26, labelLongWidth, textHeight], ...
      'ToolTipString', ...
      'Trigger for switching frame centers in MADS-PB (% of f(x_k))');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.other{4}, ...
      'Position',     [labelLeft, .23, labelLongWidth, textHeight], ...
      'ToolTipString', ...
      ['If a discrete neighbor has objective function value within this',...
       ' value of that of the incumbent, then extended polling is', ...
       ' extended polling is performed around this neighbor point']);
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.other{5}, ...
      'Position',     [labelLeft, .20, labelLongWidth, textHeight], ...
      'ToolTipString', ...
      ['If a discrete neighbor has constraint violation function value', ...
       ' within this value of that of the incumbent, then extended' ...
       ' polling is performed around this neighbor point']);
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.other{6}, ...
      'Position',     [labelLeft, .17, labelLongWidth, textHeight], ...
      'ToolTipString', ...
      ['For stochastic optimization problems, this value is the number of ', ...
       'runs to be performed for the purpose of replication.']);

% Display of NOMADm object values (Problem Name, Parameter Settings, etc.)
gui.problem      = uicontrol(gui.fig, 'String',          '', ...
                                      'FontWeight',      'bold', ...
                                      'ForegroundColor', 'red', ...
                                      'Position', [.185,.95,.265,textHeight]);
gui.searchLabel(1) = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldShortLeft, .87, fieldLongWidth,  textHeight]);
gui.searchLabel(2) = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldShortLeft, .84, fieldLongWidth,  textHeight]);
gui.pollStrategy   = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldShortLeft, .74, fieldLongWidth,  textHeight]);
gui.pollOrder      = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldShortLeft, .71, fieldLongWidth,  textHeight]);
gui.pollCenter     = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldShortLeft, .68, fieldLongWidth,  textHeight]);
gui.Term.delta     = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .64, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.Term.nIter     = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .61, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.Term.nFunc     = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .58, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.Term.time      = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .55, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.Term.nFails    = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .52, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.delta0         = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .48, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.deltaMax       = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .45, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.meshRefine     = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .42, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.meshCoarsen    = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .39, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.tolCache       = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .36, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.hmin           = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .32, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.hmax           = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .29, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.hRho           = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .26, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.ePollXiF       = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .23, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.ePollXiH       = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .20, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');
gui.nRuns          = uicontrol(gui.fig, 'ForegroundColor','blue','Position', ...
                     [fieldLongLeft,  .17, fieldShortWidth, textHeight], ...
                     'HorizontalAlignment','center');

% Additional ToolTips for NOMADm objects
set([gui.searchLabel], 'ToolTipString', ...
    'The Search screen is accessed in the MADS Menu');
set([gui.pollStrategy;gui.pollOrder;gui.pollCenter], 'ToolTipString', ...
    'Changing these values is done in the MADS Menu');
set(gui.hmin, 'ToolTipString', ...
    'This value should be kept very small');
set(gui.hmax, 'ToolTipString', ...
    'This value must be greater than the Minimum Filter Constraint Violation');
set(gui.TermFlag.delta,'ToolTipString',...
    'MADS terminates when the poll size becomes smaller than this value');
set([gui.TermFlag.nIter;gui.TermFlag.nFunc;gui.TermFlag.time; ...
     gui.TermFlag.nFails], 'ToolTipString', ...
    'Check the box to activate the corresponding termination criterion');
set(gui.Term.delta,'ToolTipString', ...
    'Keep this positive value small, as it measures solution accuracy');
set([gui.Term.nIter;gui.Term.nFunc;gui.Term.time;gui.Term.nFails], ...
    'ToolTipString', ...
    'Setting this value to ''Inf'' is the same as unchecking the box');

% Display of 2 plot axes and Stop/Resume Run pushbuttons
gui.axesHistory = axes('Parent',gui.fig, ...
   'Position',            [.52 .605 .43 .33], ...
   'FontUnits',           'normalized', ...
   'FontSize',            .07, ...
   'Box',                 'on', ...
   'Visible',             'off', ...
   'ButtonDownFcn',       gui.func.copyPlot);
gui.axesFilter = axes('Parent',gui.fig, ...
   'Position',            [.62 .13 .33 .33], ...
   'FontUnits',           'normalized', ...
   'FontSize',            .07, ...
   'Box',                 'on', ...
   'Visible',             'off', ...
   'ButtonDownFcn',       gui.func.copyPlot);
gui.stopRun = uicontrol(gui.fig, ...
   'Style',               'pushbutton', ...
   'String',              'Stop Run', ...
   'Tag',                 'StopRun', ...
   'FontWeight',          'bold', ...
   'HorizontalAlignment', 'center', ...
   'Position',            [.08 .05 .11 .045], ...
   'FontSize',            .5, ...
   'BackgroundColor',     [.8 .8 .8], ...
   'Visible',             'off', ...
   'Interruptible',       'off', ...
   'UserData',            0, ...
   'Callback',            'set(gcbo,''UserData'',1);');
gui.resumeRun = uicontrol(gui.fig, ...
   'Style',               'pushbutton', ...
   'String',              'Resume Run', ...
   'FontWeight',          'bold', ...
   'HorizontalAlignment', 'center', ...
   'Position',            [.22 .05 .11 .045], ...
   'FontSize',            .5, ...
   'BackgroundColor',     [.8 .8 .8], ...
   'Visible',             'off', ...
   'Interruptible',       'off', ...
   'Callback',            {gui.func.runMADS, 2});

% Set up Main Menu
gui.Menu.Problem = uimenu(gui.fig,'Position',1,'Label','Problem');
gui.Menu.MADS    = uimenu(gui.fig,'Position',2,'Label','MADS');
gui.Menu.Options = uimenu(gui.fig,'Position',3,'Label','Options');
gui.Menu.Session = uimenu(gui.fig,'Position',4,'Label','Session');
gui.Menu.Cache   = uimenu(gui.fig,'Position',5,'Label','Cache');
gui.Menu.Run     = uimenu(gui.fig,'Position',6,'Label','Run');
gui.Menu.Results = uimenu(gui.fig,'Position',7,'Label','Results');
gui.Menu.Help    = uimenu(gui.fig,'Position',8,'Label','Help');

% Define Problem Menu
gui.ProblemMenu.new     = uimenu(gui.Menu.Problem, 'Position',  1, ...
   'Label',       'New Optimization Problem', ...
   'Accelerator', 'N',  ...
   'Enable',      'on', ...
   'Callback',    gui.func.selectProblem);
for k = 1:length(gui_var.Types.file)
    gui.ProblemMenu.edit(k) = uimenu(gui.Menu.Problem, 'Position', k+1, ...
      'Label',    ['Edit ', gui_var.Labels.file{k}, ' File'], ...
      'Enable',   'off', ...
      'CallBack',  gui.func.editFile);
end
gui.ProblemMenu.quit = uimenu(gui.Menu.Problem, ...
   'Position',    length(gui_var.Types.file)+2, ...
   'Label',       'Quit NOMADm', ...
   'Accelerator', 'Q', ...
   'Callback',    {gui.func.quitGUI, gui_var});
set([gui.ProblemMenu.edit(1);gui.ProblemMenu.quit],'Separator','on');

% Define Options Menu
gui.OptionsMenu.scaleMenu         = uimenu(gui.Menu.Options, 'Position', 1, ...
   'Label',       'Scaling of Mesh Directions');
gui.OptionsMenu.filterMenu        = uimenu(gui.Menu.Options, 'Position', 2, ...
   'Label',       'Filter for Nonlinear Constraints');
gui.OptionsMenu.multiFiMenu       = uimenu(gui.Menu.Options, 'Position', 3, ...
   'Label',       'Run as Multifidelity Problem');
gui.OptionsMenu.sensorMenu        = uimenu(gui.Menu.Options, 'Position', 4, ...
   'Label',       'Run as Sensor Placement Problem');
gui.OptionsMenu.degeneracyMenu    = uimenu(gui.Menu.Options, 'Position', 5, ...
   'Label',       'Scheme for Handling Degenerate Linear Constraints');
gui.OptionsMenu.removeRedundancy  = uimenu(gui.Menu.Options, 'Position', 6, ...
   'Label',       'Discard Redundant Linear Constraints', ...
   'Callback',    gui.func.uitoggle);
gui.OptionsMenu.runStochastic     = uimenu(gui.Menu.Options, 'Position', 7, ...
   'Label',       'Run as Stochastic Optimization Problem', ...
   'Callback',    gui.func.uitoggle);
gui.OptionsMenu.fixCategorical    = uimenu(gui.Menu.Options, 'Position', 8, ...
   'Label',       'Fix Categorical Variables as Constant', ...
   'Callback',    gui.func.uitoggle);
gui.OptionsMenu.accelerate        = uimenu(gui.Menu.Options, 'Position', 9, ...
   'Label',       'Accelerate Convergence', ...
   'Callback',    gui.func.uitoggle);
gui.OptionsMenu.TermFlag.relative = uimenu(gui.Menu.Options, 'Position',10, ...
   'Label',       'Use Relative Termination Tolerance', ...
   'Callback',    gui.func.uitoggle);
gui.OptionsMenu.debugMenu         = uimenu(gui.Menu.Options, 'Position',11, ...
   'Label',       'Print Debugging Messages to Screen', 'Separator', 'on');
gui.OptionsMenu.saveHistory       = uimenu(gui.Menu.Options, 'Position',12, ...
   'Label',       'Save History to Text File', ...
   'Callback',    gui.func.uitoggle);
gui.OptionsMenu.plotHistoryMenu   = uimenu(gui.Menu.Options, 'Position',13, ...
   'Label',       'Plot History');
gui.OptionsMenu.plotFilter        = uimenu(gui.Menu.Options, 'Position',14, ...
   'Label',       'Plot Filter (Real-time)', ...
   'Callback',    gui.func.uitoggle);

% Define Scale SubMenu of the Options Menu choice Scale
for k = 1:length(gui_var.Labels.scale)
   gui.scaleMenu(k) = uimenu(gui.OptionsMenu.scaleMenu, ...
      'Position', k, ...
      'Label',    gui_var.Labels.scale{k}, ...
      'Callback', ['set(get(get(gcbo,''Parent''),''Children''),', ...
                   ' ''Checked'',''off''); set(gcbo,''Checked'',''on'');']);
end

% Define Filter SubMenu of the Options Menu choice Filter
for k = 1:length(gui_var.Labels.filter)
   gui.filterMenu(k) = uimenu(gui.OptionsMenu.filterMenu, ...
      'Position', k, ...
      'Label',    gui_var.Labels.filter{k}, ...
      'Callback', ['set(get(get(gcbo,''Parent''),''Children''),', ...
                   ' ''Checked'',''off''); set(gcbo,''Checked'',''on'');']);
end

% Define Filter SubMenu of the Options Menu choice Filter
for k = 1:length(gui_var.Labels.multiFi)
   gui.multiFiMenu(k) = uimenu(gui.OptionsMenu.multiFiMenu, ...
      'Position', k, ...
      'Label',    gui_var.Labels.multiFi{k}, ...
      'Callback', ['set(get(get(gcbo,''Parent''),''Children''),', ...
                   ' ''Checked'',''off''); set(gcbo,''Checked'',''on'');']);
end

% Define Sensor SubMenu of the Options Menu choice Sensor
for k = 1:length(gui_var.Labels.sensor)
   gui.sensorMenu(k) = uimenu(gui.OptionsMenu.sensorMenu, ...
      'Position', k, ...
      'Label',    gui_var.Labels.sensor{k}, ...
      'Callback', ['set(get(get(gcbo,''Parent''),''Children''),', ...
                   ' ''Checked'',''off''); set(gcbo,''Checked'',''on'');']);
end

% Define Degeneracy SubMenu of the Options Menu choice Degeneracy
for k = 1:length(gui_var.Labels.degeneracyScheme)
   gui.degeneracyMenu(k) = uimenu(gui.OptionsMenu.degeneracyMenu, ...
      'Position', k, ...
      'Label',    gui_var.Labels.degeneracyScheme{k}, ...
      'Callback', ['set(get(get(gcbo,''Parent''),''Children''),', ...
                   ' ''Checked'',''off''); set(gcbo,''Checked'',''on'');']);
end

% Define Debug SubMenu of the Options Menu choice Debug
for k = 1:length(gui_var.Labels.debug)
   gui.debugMenu(k) = uimenu(gui.OptionsMenu.debugMenu, ...
      'Position', k, ...
      'Label',    gui_var.Labels.debug{k}, ...
      'Callback', ['set(get(get(gcbo,''Parent''),''Children''),', ...
                   ' ''Checked'',''off''); set(gcbo,''Checked'',''on'');']);
end

% Define Plot History SubMenu of the Options Menu choice PlotHistory
for k = 1:length(gui_var.Labels.plotHistory)
   gui.plotHistoryMenu(k) = uimenu(gui.OptionsMenu.plotHistoryMenu, ...
      'Position', k, ...
      'Label',    gui_var.Labels.plotHistory{k}, ...
      'Callback', ['set(get(get(gcbo,''Parent''),''Children''),', ...
                   ' ''Checked'',''off''); set(gcbo,''Checked'',''on'');']);
end

% Define MADS Menu
gui.MADSMenu.editSearch     = uimenu(gui.Menu.MADS, 'Position', 1, ...
   'Label',       'Select Search Strategies', ...
   'Callback',    gui.func.search_gui);
gui.MADSMenu.editPoll       = uimenu(gui.Menu.MADS, 'Position', 2, ...
   'Label',       'Select Poll Parameters', ...
   'Callback',    gui.func.poll_gui);
gui.MADSMenu.editTermParam  = uimenu(gui.Menu.MADS, 'Position', 3, ...
   'Separator',   'on', ...
   'Label',       'Edit Termination Parameters', ...
   'UserData',    {gui_var.Labels.Parameters.term, ...
                  [gui.Term.delta;gui.Term.nIter;gui.Term.nFunc; ...
                  gui.Term.time;gui.Term.nFails]}, ...
  'Callback',     gui.func.edit_gui);
gui.MADSMenu.editMeshParam  = uimenu(gui.Menu.MADS, 'Position', 4, ...
   'Label',       'Edit Mesh Parameters', ...
   'UserData',    {gui_var.Labels.Parameters.mesh, ...
                  [gui.delta0;gui.deltaMax;gui.meshRefine; ...
                  gui.meshCoarsen;gui.tolCache]}, ...
   'Callback',    gui.func.edit_gui);
gui.MADSMenu.editOtherParam = uimenu(gui.Menu.MADS, 'Position', 5, ...
   'Label',       'Edit Filter/MVP/Stochastic Parameters', ...
   'UserData',    {gui_var.Labels.Parameters.other, ...
                  [gui.hmin; gui.hmax; gui.hRho; ...
                   gui.ePollXiF; gui.ePollXiH; gui.nRuns]}, ...
   'Callback',    gui.func.edit_gui);

% Define Session Menu
gui.SessionMenu.load   = uimenu(gui.Menu.Session, 'Position', 1, ...
   'Label',       'Load Options from Session File', ...
   'Enable',      'off', ...
   'Callback',    gui.func.loadSession);
gui.SessionMenu.save   = uimenu(gui.Menu.Session, 'Position', 2, ...
   'Label',       'Save Options to Session File', ...
   'Enable',      'off', ...
   'Callback',    gui.func.saveSession);
gui.SessionMenu.reset  = uimenu(gui.Menu.Session, 'Position', 3, ...
   'Label',       'Reset Options to Defaults', ...
   'Callback',    gui.func.resetSession);
gui.SessionMenu.delete = uimenu(gui.Menu.Session, 'Position', 4, ...
   'Label',       'Delete Session File', ...
   'Enable',      'off', ...
   'UserData',   {'Session File',gui_var.FileExt.S,gui.SessionMenu.load},...
   'Callback',    gui.func.deleteFile);

% Define Cache Menu
gui.CacheMenu.load   = uimenu(gui.Menu.Cache, 'Position', 1, ...
   'Label',       'Use Pre-Existing Cache File', ...
   'Enable',      'off', ...
   'Callback',    gui.func.uitoggle);
gui.CacheMenu.count  = uimenu(gui.Menu.Cache, 'Position', 2, ...
   'Label',       'Count Cache Points as Function Calls', ...
   'Enable',      'off', ...
   'Callback',    gui.func.uitoggle);
gui.CacheMenu.save   = uimenu(gui.Menu.Cache, 'Position', 3, ...
   'Label',       'Save Run Results to Cache File', ...
   'Accelerator', 'S', ...
   'Enable',      'off', ...
   'Callback',    gui.func.saveCache);
gui.CacheMenu.delete = uimenu(gui.Menu.Cache, 'Position', 4, ...
   'Label',       'Delete Pre-Existing Cache File', ...
   'Enable',      'off', ...
   'UserData',    {'Cache File', gui_var.FileExt.C, ...
                  [gui.CacheMenu.load; gui.CacheMenu.count]}, ...
   'Callback',    gui.func.deleteFile);
gui.CacheMenu.countNotInX = uimenu(gui.Menu.Cache, 'Position', 5, ...
   'Label',       'Add X-infeasible points to Cache', ...
   'Visible',     'off', ...
   'Callback',    gui.func.uitoggle);

% Define Run Menu
gui.RunMenu.exec         = uimenu(gui.Menu.Run,'Position',1, ...
   'Label',       'Execute Next Run', ...
   'Accelerator', 'R', ...
   'Enable',      'off', ...
   'Callback',    {gui.func.runMADS, 0});
gui.RunMenu.resume       = uimenu(gui.Menu.Run,'Position',2, ...
   'Label',       'Resume a Stopped Run', ...
   'Enable',      'off', ...
   'Callback',    {gui.func.runMADS, 2});
gui.RunMenu.restart      = uimenu(gui.Menu.Run,'Position',3, ...
   'Label',       'Restart Run from Current Point', ...
   'Enable',      'off', ...
   'Callback',    {gui.func.runMADS, 1});
gui.RunMenu.oneIteration = uimenu(gui.Menu.Run,'Position',4, ...
   'Separator',   'on', ...
   'Label',       'Run One Iteration at a Time', ...
   'Checked',     'off', ...
   'Enable',      'off', ...
   'Callback',    gui.func.uitoggle);
gui.RunMenu.execFeasible = uimenu(gui.Menu.Run,'Position',5, ...
   'Label',       'Run Only Until Feasible', ...
   'Checked',     'off', ...
   'Enable',      'off', ...
   'Callback',    gui.func.uitoggle);
gui.RunMenu.clear1       = uimenu(gui.Menu.Run,'Position',6, ...
   'Separator',   'on', ...
   'Label',       'Reset Previous Run', ...
   'Enable',      'off', ...
   'Callback',    {gui.func.clearRuns,0});
gui.RunMenu.clear        = uimenu(gui.Menu.Run,'Position',7, ...
   'Label',       'Clear Previous Runs', ...
   'Enable',      'off', ...
   'Callback',    {gui.func.clearRuns,1});

% Define Results Menu
set(gui.Menu.Results,'Visible','off');
for k = 1:gui_var.runMax
   gui.ResultsMenu(k) = uimenu(gui.Menu.Results,'Position',k, ...
      'Label',          ['View Run # ' int2str(k)], ...
      'Visible',        'off', ...
      'ForegroundColor',gui_var.Types.plotColors(k), ...
      'Accelerator',    int2str(mod(k,gui_var.runMax)), ...
      'Callback',       gui.func.results_gui);
end

% Define Help Menu
uimenu(gui.Menu.Help,'Position',1, ...
   'Label',      'NOMADm Help', ...
   'Accelerator','H', ...
   'Visible',    onoff{2-~~exist(gui_var.HelpDoc.nomadm,'file')}, ...
   'UserData',   gui_var.HelpDoc.nomadm,...
   'Callback',   gui.func.getHelp);
uimenu(gui.Menu.Help,'Position',2, ...
   'Label',      'DACE Toolbox Help', ...
   'Visible',    onoff{2-~~exist(gui_var.HelpDoc.dace,'file')}, ...
   'UserData',   gui_var.HelpDoc.dace,...
   'Callback',   gui.func.getHelp);
uimenu(gui.Menu.Help,'Position',3, ...
   'Label',      'N-W Toolbox Help', ...
   'Visible',    onoff{2-~~exist(gui_var.HelpDoc.nw,'file')}, ...
   'UserData',   gui_var.HelpDoc.nw,...
   'Callback',   gui.func.getHelp);
uimenu(gui.Menu.Help,'Position',4, ...
   'Label',      'RBF Toolbox Help', ...
   'Visible',    onoff{2-~~exist(gui_var.HelpDoc.rbf,'file')}, ...
   'UserData',   gui_var.HelpDoc.rbf,...
   'Callback',   gui.func.getHelp);
uimenu(gui.Menu.Help,'Position',5, ...
   'Separator',  'on', ...
   'Label',      'View List of Version Changes', ...
   'Visible',    onoff{2-~~exist(gui_var.HelpDoc.changes,'file')}, ...
   'UserData',   gui_var.HelpDoc.changes,...
   'Callback',   gui.func.getHelp);
uimenu(gui.Menu.Help,'Position',6, ...
   'Label',      'View GNU Public License', ...
   'Visible',    onoff{2-~~exist(gui_var.HelpDoc.license,'file')}, ...
   'UserData',   gui_var.HelpDoc.license,...
   'Callback',   gui.func.getHelp);
uimenu(gui.Menu.Help,'Position',7, ...
   'Separator',  'on', ...
   'Label',      'About NOMADm', ...
   'Callback',   {gui.func.getAbout,gui_var.Labels.version,...
                  gui_var.Labels.coolCats});
end

%*******************************************************************************
% BEGINNING OF SIMPLE NOMADm GUI CALLBACK FUNCTIONS
%*******************************************************************************
%*******************************************************************************
% selectProblem: Choose an optimization problem.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Problem-->New Optimization Problem)
% Calls:     clearRuns, resetSession
% VARIABLES:
%  h                 = handle of calling object
%  pFile             = filename of current optimization problem
%  pPath             = pathname of current optimization problem
%  gui_var           = structure of all GUI variables
%    .problemExt     =   filename extension of optimization problem 
%    .problem        =   name of current optimization problem
%    .FileExt        =   filename suffixes
%      .O            =     Omega file suffix
%      .I            =     initial points file suffix
%      .N            =     discrete neighbors file suffix
%      .P            =     user parameter file suffix
%      .C            =     Cache file suffix
%      .S            =     Session file suffix
%    .noGrad         =   flag indicating initial availability of gradients
%  gui               = structure of all GUI object handles
%    .runStatus      =   handle for text bar at bottom of figure window
%    .RunMenu        =   handles for Run menu items
%    .ProblemMenu    =   handles for Problem menu items
%    .MADSMenu       =   handles for MADS menu items
%    .SessionMenu    =   handles for Session menu items
%    .CacheMenu      =   handles for CacheMenu items
%  nameProblem       = name of current optimization problem
%  isConstrained     = flags problem as having linear constraints
%  hasInitGuess      = flags problem as having an initial guess
%  hasUserParam      = flags problem as having user parameters
%  isMVP             = flags problem as being an MVP
%  existCache        = flags problem as having an existing Cache file
%  existSession      = flags problem as having an existing Session file
%  pLabels           = strings of Poll strategy labels
%*******************************************************************************
function selectProblem(h,~)  %#ok
[pFile,pPath] = uigetfile({'*.m',        'Matlab M-files (*.m)';      ...
                           '*.f; *.for', 'Fortran files (*.f,*.for)'; ...
                           '*.c; *.C',   'C/C++ files (*.c,*.C)'},    ...
                           'Select an Optimization Problem file');
if ~pFile, return, end

% Reset all variables, retrieve GUI application data, and group menu handles
clearRuns(h,[],1);
resetSession(h,[]);
gui_var      = getappdata(gcbf,'gui_var');
gui          = getappdata(gcbf,'gui');
Opts.Run     = [gui.RunMenu.exec; gui.RunMenu.oneIteration];
Opts.Run     = [Opts.Run; gui.RunMenu.execFeasible];
Opts.Cache   = [gui.CacheMenu.load; gui.CacheMenu.count; gui.CacheMenu.delete];
Opts.Session = [gui.SessionMenu.load; gui.SessionMenu.delete];
if ~isdeployed, cd(pPath); end

% Load New Problem
[nameProblem, gui_var.problemExt] = strtok(pFile,'.');
set(gui.problem,   'String',nameProblem);
set(gui.runStatus, 'String','No runs performed');

% Allow editing only of files that exist
fld   = fieldnames(gui_var.FileExt);
onoff = {'on','off'};
on{1} =  'on';
for k = 2:length(fld)
   on{k} = onoff{1+~exist([nameProblem,gui_var.FileExt.(fld{k})],'file')}; %#ok
end
for k = 1:length(gui.ProblemMenu.edit)
   set(gui.ProblemMenu.edit(k), 'Enable', on{k});
end

% Enable/disable grouped menu items, as appropriate
set(gui.CacheMenu.countNotInX, 'Visible', on{3});
set(Opts.Cache,                'Enable',  on{7});
set(Opts.Session,              'Enable',  on{8});
set(Opts.Run,                  'Enable', 'on');
set(gui.SessionMenu.save,      'Enable', 'on');

% Allow gradient options if functions file returns enough arguments
try
   gui_var.noGrad = (abs(nargout(nameProblem)) <= 2);
catch
   gui_var.noGrad = 0;
end

setappdata(gcbf,'gui_var',gui_var);
end

%*******************************************************************************
% editFile: Edit an optimization problem file.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Problem-->Edit XXXXX File)
% VARIABLES:
%  h             = handle of calling object
%  k             = Problem Menu position number of selected file
%  ext           = filename extension of file to be edited
%  filename      = full name of file to be edited
%  gui_var       = structure of all GUI variables
%    .problemExt =   filename extension of optimization problem 
%    .Types.file =   strings of file suffixes
%  gui.problem   = name of current optimization problem
%*******************************************************************************
function editFile(~,~) %#ok

[h,fig] = gcbo;
gui_var = getappdata(fig,'gui_var');
gui     = getappdata(fig,'gui');

k = get(h,'Position') - 1;
if (k == 1)
   ext = gui_var.problemExt;
else
   ext = '.m';
end
filename = [get(gui.problem,'String'),gui_var.Types.file{k},ext];
edit(filename);
clear(filename);
end

%*******************************************************************************
% quitGUI: End program session.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Problem-->Quit NOMADm)
%            nomadm_gui (Close Request Function for figure window)
% VARIABLES:
%  gui_var = structure of all NOMADm variables
%    .path = path from which NOMADm was launched
%*******************************************************************************
function quitGUI(~,~, gui_var) %#ok

if ~isdeployed, cd(gui_var.path); end
if ishandle(gcbf), delete(gcbf); end
end

%*******************************************************************************
% selectPoll: Select a Poll Strategy, Poll Center, or Poll Order.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: MADS-->Select Poll Directions), 
%            nomadm_gui (callback: MADS-->Select Poll Center),
%            nomadm_gui (callback: MADS-->Select Poll Order)
% VARIABLES:
%  h              = handle of calling object
%  type           = string containing "type", "center", or "order"
%  gui_var.Choice = user choices
%  gui            = structure of NOMADm object handles
%*******************************************************************************
function selectPoll(h,~,type) %#ok

% Update Poll option and screen
set(get(get(h,'Parent'),'Children'),'Checked','off');
set(h,'Checked','on');
gui_var = getappdata(gcbf,'gui_var');
gui_var.Choice.Poll.(type) = get(h,'Position');
set(get(h,'UserData'),'String',get(h,'Label'));
setappdata(gcbf,'gui_var',gui_var);
end

%*******************************************************************************
% uitoggle: Toggle checked status of uimenu object.
% ------------------------------------------------------------------------------
% Called by: many functions
% VARIABLES:
%   handle = handle of uimenu object
%   state  = status of "checked property (1=on, 2=off)
%*******************************************************************************
function uitoggle(h,~) %#ok
offon = {'off','on'};
set(h,'checked',offon{~strcmp(get(h,'checked'),'on')+1});
end

%*******************************************************************************
% loadSession: Load session options from previously saved file.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Session-->Load Options from Session File)
% Calls:     updateScreen
% VARIABLES:
%  sessionFile     = name of Session file
%  gui_var         = structure of NOMADm variables
%    .FileExt.S    = filename extension of session file
%    .Choice       = current user choices
%    .Options      = current option settings
%  gui.problem     = object handle for optimization problem name
%  Session.Choice  = previously saved user choices 
%  Session.Options = previously saved option settings
%*******************************************************************************
function loadSession(~,~) %#ok

gui_var = getappdata(gcf,'gui_var');
gui     = getappdata(gcf,'gui');

sessionFile = [get(gui.problem,'String'), gui_var.FileExt.S];
if exist(sessionFile,'file')
   load(sessionFile,'Session');
   gui_var.Choice  = Session.Choice;
   gui_var.Options = Session.Options;
   updateScreen(gui_var.Choice,gui_var.Options);
   setappdata(gcf,'gui_var',gui_var);
end
end

%*******************************************************************************
% resetSession: Reset MADS Parameters to Default values.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Session-->Reset Options to Defaults),
%            nomadm, selectProblem
% Calls:     updateScreen
% VARIABLES:
%  gui_var.Choice           = current user choices
%  gui_var.Options          = current option settings
%  gui_var.Defaults.Choice  = default user choice settings
%  gui_var.Defaults.Options = default options settings
%*******************************************************************************
function resetSession(~,~)

gui_var = getappdata(gcf,'gui_var');
gui_var.Choice  = gui_var.Defaults.Choice;
gui_var.Options = gui_var.Defaults.Options;

updateScreen(gui_var.Choice,gui_var.Options);
setappdata(gcf,'gui_var',gui_var);
end

%*******************************************************************************
% saveSession: Save session options to file for future retrieval.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Session-->Save Options to Session File)
% Calls:     loadMADS
% VARIABLES:
%  sessionFile            = name of Session file
%  gui_var                = structure of all GUI variables
%    .FileExt.S           =   filename extension of session file
%    .Choice              =   current user choices
%    .Options             =   current option settings
%  gui                    = structure of GUI object handles
%    .problem             =   name of current optimization problem
%    .SessionMenu         =   handles for Session menu items
%    .OptionsMenu         =   handles for Options menu items
%    .runStatus           =   handle for figure window status bar
%  Session                = previously saved options and parameters
%    .Problem             =   optimization problem data
%    .Choice              =   user choices 
%    .Options             =   option settings
%      .Term              =     termination criteria
%        .iter            =       number of iterations
%        .func            =       number of function evaluations
%        .time            =       CPU time
%        .fails           =       number of consecutive Poll failures
%      .TermFlag.relative =     termination relative to initial mesh size
%*******************************************************************************
function saveSession(~,~) %#ok

gui     = getappdata(gcbf,'gui');
gui_var = getappdata(gcbf,'gui_var');

sessionFile = [get(gui.problem,'String'), gui_var.FileExt.S];
[Session.Problem,Session.Options] = loadMADS(gui,gui_var);

field = fieldnames(gui.Term);
for k = 1:length(fieldnames(gui.Term))
   Session.Options.Term.(field{k}) = get(gui.Term.(field{k}), 'String');
end
Session.Options.TermFlag.relative = ...
                strcmp(get(gui.OptionsMenu.TermFlag.relative,'Checked'),'on');
Session.Choice = gui_var.Choice;
if length(Session.Options.Search) <= Session.Options.nSearches
   Session.Options.Search = [Session.Options.Search, ...
                     gui_var.Options.Search(gui_var.Options.nSearches+1:end)];
end
save(sessionFile,'Session');
set(gui.runStatus,'String', ...
   ['Session Options saved to file, ',sessionFile]);
set([gui.SessionMenu.load; gui.SessionMenu.delete],'Enable','on');
end

%*******************************************************************************
% saveCache: Save MADS run results in a .mat file for later retrieval.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Cache-->Save Run Results to Cache File),
%            loadMADS
% VARIABLES:
%  gui             = structure of all GUI object handles
%    .fig          =   handle for the GUI figure window
%    .problem      =   name of current optimization problem
%    .runStatus    =   handle for the GUI figure window status bar
%    .CacheMenu    =   handles for Cache menu items
%  gui_var         = structure of all GUI variables
%    .runCount     =   MADS run number
%    .FileExt.C    =   default Cache filename suffix 
%  CName           = name of Cache file
%  RunSet(1).Cache = structure of MADS run data
%  Cache           = storage of Cache for use by MADS
%*******************************************************************************
function saveCache(~,~)

if isappdata(gcbf,'RunSet')
   gui     = getappdata(gcbf,'gui');
   gui_var = getappdata(gcbf,'gui_var');
   RunSet  = getappdata(gcbf,'RunSet');
   setptr(gcf,'watch');
   set(gui.runStatus, 'String', ...
       ['Run # ',int2str(gui_var.runCount),' Cache being saved']);
   CName = [get(gui.problem, 'String'), gui_var.FileExt.C];
   Cache = RunSet(1).Cache; %#ok
   save(CName,'Cache');
   set([gui.CacheMenu.load;gui.CacheMenu.count;gui.CacheMenu.delete], ...
       'Enable','on');
   set(gui.runStatus, 'String', ...
       ['Run # ',int2str(gui_var.runCount),' Cache saved in ',CName]);
   setptr(gcbf,'arrow');
else
   error('nomadm:saveCache','Cannot save.  Cache not found');
end
end

%*******************************************************************************
% deleteFile: Delete a Session or Cache file.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Session-->Delete Session File), 
%            nomadm_gui (callback: Cache-->Delete Pre-Existing Cache File)
% VARIABLES:
%  h              = handle of calling object
%  param          = user data that identifies file to be deleted
%  fileID         = string identifing file to be deleted
%  fileExt        = filename suffix of file to be deleted
%  handles        = handles of menu choices to be disabled
%  file           = name of file to be deleted
%  deleteFile     = logical for deleting file
%  gui            = structure of all GUI object handles
%    .runStatus   =   object handle for figure window status bar
%    .problem     =   name of current optimization problem
%*******************************************************************************
function deleteFile(h,~) %#ok

gui     = getappdata(gcbf,'gui');
param   = get(h,'UserData');
[fileID,fileExt,handles] = deal(param{:});
file    = [get(gui.problem,'String'),fileExt];
deleteFile = questdlg(['Are you sure you want to delete ',file,'?'], ...
                      ['Delete ',fileID,'?'],'No');
if strcmp(deleteFile, 'Yes')
   delete(file);
   set(gui.runStatus,'String', [fileID,', ',file,', has been deleted']);
   set([handles; h],'Enable','off');
end
end

%*******************************************************************************
% clearRuns: Clear figure window and reset all the appropriate variables.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Run-->Clear Previous Runs), selectProblem
% VARIABLES:
%  gui               = structure of all GUI object handles
%    .fig            =   handle for GUI figure window
%    .stopRun        =   handle for Stop Run pushbutton
%    .resumeRun      =   handle for Resume Run pushbutton
%    .axesHistory    =   handle for History plot
%    .axesFilter     =   handle for Filter plot
%    .ResultsMenu    =   handles for Results menu items
%    .RunMenu        =   handles for Run menu items
%    .CacheMenu      =   handles for Cache menu items
%    .runStatus      =   handle for figure window status bar
%  gui_var           = structure of all GUI variables
%    .runMax         =   maximum nmber of MADS runs
%    .runCount       =   MADS run counter
%  RunSet            = structure of MADS run data
%  clearall          = flag for clearing all runs or resetting only previous one
%*******************************************************************************
function clearRuns(~,~,clearall)

gui     = getappdata(gcbf,'gui');
gui_var = getappdata(gcbf,'gui_var');

% Test for clearing all or only one run
if clearall || ~gui_var.runCount

   % Delete application data
   if isappdata(gcbf,'RunSet'), rmappdata(gcbf,'RunSet'); end

   % Reset properties of gui object handles
   hEnableOff  = [gui.RunMenu.clear; gui.RunMenu.clear1; gui.CacheMenu.save];
   hVisibleOff = [gui.Menu.Results;  gui.ResultsMenu'; gui.axesHistory; ...
                  get(gui.axesHistory,'Children')];
   set(hEnableOff, 'Enable', 'off');
   set(hVisibleOff,'Visible','off');
   set(gui.axesHistory,'NextPlot','replace');
   set(gui.runStatus,  'String',  'All runs cleared');

   % Reset the run counter and reset everything else
   gui_var.runCount = 0;
   setappdata(gcbf,'gui_var',gui_var);
   cla; clc;
else
   set(gui.runStatus, 'String', 'Previous run cleared');
end

% Actions common to "Clear all runs" or "Reset previous run" options
if isappdata(0,'CACHE'),     rmappdata(0,'CACHE');     end
if isappdata(0,'sCACHE'),    rmappdata(0,'sCACHE');    end
if isappdata(0,'surrogate'), rmappdata(0,'surrogate'); end
set([gui.axesFilter;    get(gui.axesFilter, 'Children')],'Visible','off');
set([gui.RunMenu.resume;gui.RunMenu.restart],            'Enable', 'off');      
set([gui.stopRun;       gui.resumeRun],                  'Visible','off');      
setptr(gcbf,'arrow');
end

%*******************************************************************************
% getHelp: View a Help file.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Help-->NOMADm Help),
%            nomadm_gui (callback: Help-->DACE Toolbox Help),
%            nomadm_gui (callback: Help-->NW Toolbox Help),
%            nomadm_gui (callback: Help-->CMA-ES Toolbox Help),
%            nomadm_gui (callback: Help-->View List of Version CHanges),
%            nomadm_gui (callback: Help-->View GNU Public License),
% VARIABLES:
%  h         = handle of calling object
%  helpfile  = name of the help file to be displayed
%  errormsg  = message to display when error is flagged
%  errorflag = error flag for viewing failure
%*******************************************************************************
function getHelp(h,~) %#ok

helpfile = which(get(h,'UserData'));
if isempty(helpfile)
   errormsg = ['Error: Help file not found on the Matlab path.  ', ...
               'Use "File-->Set Path" to add location to Matlab path.'];
   errordlg(errormsg);
end
errorflag = web(['file:///', helpfile], '-browser');
if (errorflag), errordlg('Error: Browser or help file not found'); end
end

%*******************************************************************************
% getAbout: View author information.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Help-->About NOMADm)
% VARIABLES:
%    h            = handle of calling object
%    pversion     = current version number of this software
%    contributors = list of people who have contributed to this software
%    today        = string containing today's date
%*******************************************************************************
function getAbout(h,~,pversion,contributors) %#ok
today = date;
msgbox([{['NOMADm, version ',pversion]}, ...
        {['Copyright (c) 2001-',today(end-3:end),' by Mark A. Abramson']}, ...
        {''},{'Special Thanks to:'},{''},contributors], get(h,'Label'));
end

%*******************************************************************************
% copyPlot: Copy history or filter plot to a new screen.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (ButtonDown Function: History and Filter Plots)
% VARIABLES:
%  h  = handle of calling object
%  h1 = handle to new figure window
%  h2 = handle to new axes on figure window
%*******************************************************************************
function copyPlot(~,~) %#ok

h1 = figure;
h2 = copyobj(gcbo,h1);
set(h2,'Position','default','ButtonDownFcn','');
end

%*******************************************************************************
% END OF SIMPLE NOMADm GUI CALLBACK FUNCTIONS
%*******************************************************************************
%*******************************************************************************
% runMADS:  Assign GUI input fields to MADS variables and run MADS.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Run-->Execute Next Run)
%            nomadm_gui (callback: Run-->Resume a Stopped Run)
%            nomadm_gui (callback: Run-->Restart Run from Current Point)
% Calls:     loadMADS, saveCache, mads
% VARIABLES:
%  restart           = flag for resuming previous run
%  gui_var           = structure containing all GUI variables
%    .runCount       =   current MADS run number
%    .runMax         =   maximum allowed MADS run number
%    .Options.runUntilFeasible = flag for running until feasible point is found
%  gui.func          = structure of all NOMADm function handles
%  gui               = structure of GUI object handles
%    .fig            =   handle for the main figure window
%    .stopRun        =   handle for Stop Run pushbutton
%    .resumeRun      =   handle for Resume Run pushbutton
%    .axesFilter     =   handle for Filter plot
%    .axesHistory    =   handle for History plot
%    .RunMenu        =   handles for Run menu items
%    .ResultsMenu    =   handles for Results menu items
%    .CacheMenu      =   handles for Cache menu items
%    .runStatus      =   handle for figure window status bar
%  RunSet            = global variable used to store MADS output data
%    .BestF          =   best feasible solution found by MADS
%    .BestI          =   best infeasible solution found by MADS
%    .RunStats       =   structure of MADS run statistics
%      .delta        =     current poll size
%      .time         =     current CPU time expended
%      .Cache        =   structure of all processed MADS iterates
%  Problem           = structure containing optimization problem data
%    .File.I         =   name of initial points file
%  Options           = structure containing MADS parameters
%    .loadCache      =   flag for loading a pre-existing Cache file
%    .countCache     =   flag for counting Cache points as function calls
%    .delta0         =   initial poll size
%    .runCount       =   current MADS run number
%  iterate0          = initial iterate
%  Param             = structure of output from user-provided Parameter file
%  BestF             = best feasible iterate found by MADS
%  BestI             = least infeasible iterate found by MADS
%  RunStats          = structure containing MADS run statistics
%*******************************************************************************
function runMADS(~,~,restart) %#ok

% Get application data
fig     = findobj(0,'Tag','NOMADm_GUI');
gui_var = getappdata(fig,'gui_var');
gui     = getappdata(fig,'gui');
if isappdata(fig,'RunSet')
   RunSet = getappdata(fig,'RunSet');
end

% Set flag for running only until a feasible point is found
if (restart == -1)
   gui_var.Options.runUntilFeasible = 1;
end

% lasterr('NOMADm interrupted by user');
try

   % Change figure window, as appropriate
   if (gui_var.runCount >= gui_var.runMax)
      error('nomadm:runMADS','Too many MADS runs without clearing.');
   end
   setptr(fig,'watch');
   set(gui.stopRun,   'Visible','on','UserData',0);
   set(gui.resumeRun, 'Visible','off');
   set([gui.RunMenu.clear1; gui.RunMenu.clear],'Enable','on');
   set([gui.axesFilter; get(gui.axesFilter,'Children')], 'Visible', 'off');
   set(gui.axesFilter,'ButtonDownFcn','');
   drawnow;

   % Load MADS input data
   [Problem,Options] = loadMADS(gui,gui_var);
   Options.Search(Options.nSearches+1:end) = [];

   % Get optional Parameter data from user-provided file
   if (exist(Problem.File.P,'file') == 2)
      Problem.Param = feval(Problem.File.P);
      setappdata(0,'PARAM',Problem.Param);
   end

   % Get initial point
   if isfield(Problem,'Param') && isfield(Problem.Param,'iterate0')
      iterate0 = Problem.Param.iterate0;
   else
      if (restart > 0)
         iterate0 = [RunSet(gui_var.runCount).BestF, ...
                     RunSet(gui_var.runCount).BestI];
      elseif (~restart && exist(Problem.File.I,'file'))
         iterate0 = feval(Problem.File.I);
      else
         error('nomadm:runMADS',['Cannot find: ', Problem.File.I, '.m']);
      end
   end

   % Set up and run MADS algorithm and store output
   if (restart == 2)
      saveCache(1,[]);
      set(gui.axesHistory, 'NextPlot','replacechildren');
      Options.loadCache  = 1;
      Options.countCache = 1;
      Options.plotColor  = gui_var.Types.plotColors(gui_var.runCount);
      Options.delta0     = RunSet(gui_var.runCount).RunStats.delta;
      RunStats           = RunSet(gui_var.runCount).RunStats;
      set(gui.runStatus,'String',['Resuming Run # ',int2str(gui_var.runCount)]);
      Problem.File.H = [Problem.File.H, int2str(gui_var.runCount),'.txt'];
      [BestF,BestI,RunStats,RunSet(1).Cache] = mads(Problem,iterate0, ...
                                                    Options,RunStats);
      delete(Problem.File.C);
      set([gui.CacheMenu.load; gui.CacheMenu.count; gui.CacheMenu.delete], ...
          'Enable','off');
   else
      set(gui.runStatus,'String', ...
                       ['Processing Run # ',int2str(gui_var.runCount+1)]);
      Problem.File.H = [Problem.File.H, int2str(gui_var.runCount+1),'.txt'];
      [BestF,BestI,RunStats,RunSet(1).Cache] = mads(Problem,iterate0,Options);
      gui_var.runCount = gui_var.runCount+1;
   end

   % Perform any user-defined post-processing (must have argument)
   if (exist(Problem.File.P,'file') == 2) && (nargin(Problem.File.P) < 0)
      Param = feval(Problem.File.P,BestF);
      setappdata(0,'PARAM',Param);
      if isfield(Param,'BestF')
         BestF = Param.BestF;
      end
      if isfield(Param,'BestI')
         BestF = Param.BestI;
      end
   end

   % Store results
   RunSet(gui_var.runCount).BestF    = BestF;
   RunSet(gui_var.runCount).BestI    = BestI;
   RunSet(gui_var.runCount).RunStats = RunStats;
   
   % Sound train whistle at termination if run takes longer than one minute
   if (RunStats.time > 60)
      try
         load('train','y');
         sound(y)
      catch
      end
   end

% Perform these tasks if error in MADS run
catch failMADS
   set(gui.runStatus,'String',['Run # ',int2str(gui_var.runCount+1),' failed']);
   set([gui.axesFilter; get(gui.axesFilter, 'Children')],'Visible','off');
   if (gui_var.runCount == 0)
      set([gui.axesHistory; get(gui.axesHistory,'Children')],'Visible','off');
   else
      set(gui.axesHistory,'ButtonDownFcn',gui.func.copyPlot);
   end
   RunSet(1).Cache = [];
   if isappdata(0,'CACHE')
      rmappdata(0,'CACHE');
   end
   errorMsg = failMADS.getReport('basic','hyperlinks','off');
   errordlg(errorMsg,'MADS Runtime Error','modal');
   beep;
   rethrow(failMADS);
end

% Change figure window, as appropriate
set(gui.axesHistory, 'NextPlot','add');
set(gui.axesHistory, 'ButtonDownFcn',gui.func.copyPlot);
set(gui.axesFilter,  'ButtonDownFcn',gui.func.copyPlot);

if (gui_var.runCount > 0)
   set([gui.CacheMenu.save;gui.RunMenu.resume;gui.RunMenu.restart],'Enable','on');
   set([gui.Menu.Results; gui.ResultsMenu(gui_var.runCount)],'Visible','on');
   set(gui.runStatus,'String',['Run # ',int2str(gui_var.runCount),' complete']);
end
set(gui.stopRun,   'Visible','off','UserData',0);
set(gui.resumeRun, 'Visible','on');
setptr(fig,'arrow');
setappdata(fig,'gui_var',gui_var);
setappdata(fig,'RunSet',RunSet);

end

%*******************************************************************************
% loadMADS:  Assign GUI input fields to MADS variables.
% ------------------------------------------------------------------------------
% Called by: saveSession, runMADS
% Calls:     nomadm_compile
% VARIABLES:
%  Problem              = structure containing optimization problem data
%    .File              =   structure of problem file names
%      .F               =   name of functions file
%      .O               =   name of Omega file
%      .I               =   name of initial points file
%      .N               =   name of discrete neighbors file
%      .C               =   name of Cache File
%      .nameCache       =   name of the base workspace Cache variable
%      .fType           =   type of functions file (M=MATLAB,F=FORTRAN,C=C)
%  Options              = structure containing MADS parameters
%    .nSearches         =   number of Search types used
%    .Search(n)         =   structure of Search parameters
%      .type            =     string identifying the type of Search
%    .Poll.type         =   string identifying selected Poll strategy
%    .Poll.order        =   string identifying selected Poll order strategy
%    .Poll.center       =   integer identifying selected Poll center
%    .Poll.complete     =   turns on/off complete Polling
%    .EPoll.completeN   =   turns on/off complete discrete neighbor Polling
%    .EPoll.complete    =   turns on/off complete Extended Polling
%    .EPoll.fTrigger    =   f-value Extended Poll trigger
%    .EPoll.hTrigger    =   h-value Extended Poll trigger
%    .loadCache         =   flag for loading a pre-existing Cache file
%    .countCache        =   flag for counting Cache points as function calls
%    .useFilter         =   use filter for nonlinear constraints
%    .degeneracyScheme  =   scheme for handling degenerate linear constraints
%    .removeRedundancy  =   remove redundant linear constraints
%    .runStochastic     =   flag for running as a stochastic problem
%    .fixCategorical    =   flag for holding categorical variables constant
%    .accelerate        =   flag for accelerating mesh refinement
%    .scale             =   flag for scaling mesh directions
%    .saveHistory       =   flag for saving history to a text file
%    .plotHistory       =   turns on/off a history plot
%    .plotFilter        =   turns on/off a real-time filter plot
%    .runOneIteration   =   flag for running one MADS iteration at a time
%    .runUntilFeasible  =   flag for running MADS only until feasible
%    .runCount          =   MADS run counter
%    .hplothandle       =   handle for history plot axes
%    .fplothandle       =   handle for filter plot axes
%    .stophandle        =   handle for Stop Run pushbutton
%    .delta0            =   initial poll size
%    .deltaMax          =   maximum poll size
%    .meshRefine        =   mesh refinement factor
%    .meshCoarsen       =   mesh coarsening factor
%    .tolCache          =   tolerance for flagging point as being in Cache
%    .hmin              =   minimum h-value of an infeasible point
%    .hmax              =   maximum h-value of a filter point
%    .Term              =   substructure containing MADS termination criteria
%      .delta           =     poll size parameter
%      .iter            =     maximum number of MADS iterations
%      .func            =     maximum number of function evaluations
%      .time            =     maximum CPU time
%      .fails           =     maximum number of consecutive Poll failures
%    .TermFlag          =   substructure of termination criteria on/off switches
%      .iter            =     turns on/off number of iterations
%      .nFunc           =     turns on/off number of function evaluations
%      .time            =     turns on/off CPU time
%      .nFails          =     turns on/off number of consecutive Poll failures
%      .relative        =     computes termination delta relative to .delta0
%  gui_var              = structure containing all GUI variables
%    .problemExt        =   filename extension of optimization problem 
%    .Options           =   current Options values
%    .Types             = lists of possible type
%      .Search          =   list of possible Search types
%      .poll            =   list of possible Poll strategies
%      .pollOrder       =   list of possible Poll order strategies
%    .Choice            = user choices
%      .pollStrategy    =   selected Poll strategy
%      .pollCenter      =   selected Poll center
%      .pollOrder       =   selected Poll order strategy
%  gui                  = structure of GUI object handles
%    .runStatus         =   handle for GUI figure window status bar
%    .problem           =   name of current optimization problem
%    .CacheMenu         =   handles for Cache menu items
%    .MADSMenu          =   handles for MADS menu items
%    .OptionsMenu       =   handles for Options menu items
%    .RunMenu           =   handles for Run menu items
%    .scaleMenu         =   handles for scaling submenu items
%    .filterMenu        =   handles for filter submenu items
%    .degeneracyMenu    =   handles for degeneracy submenu items
%    .delta0            =   current initial poll size
%    .deltaMax          =   current maximum poll size
%    .meshRefine        =   current mesh refinement factor
%    .meshCoarsen       =   current mesh coarsening factor
%    .tolCache          =   current Cache tolerance
%    .hmin              =   current minimum infeasible h-value
%    .hmax              =   current maximum filter h-value
%    .ePollXiF          =   current f-value Extended Poll trigger
%    .ePollXiH          =   current h-value Extended Poll trigger
%    .Term              =   structure of handles for current termination criteria
%      .delta           =     current poll size
%      .nIter           =     current maximum number of MADS iterations
%      .nFunc           =     current maximum number of function evaluations
%      .time            =     current maximum CPU time
%      .nFails          =   current max number of consecutive Poll failures
%    .TermFlag          =   structure of handles for termination checkboxes
%      .nIter           =     handle for number of iterations
%      .nFunc           =     handle for number of function evaluations
%      .time            =     handle for CPU time
%      .nFails          =     handle for number of consecutive Poll failures
%  nameProblem          = name of the optimization problem to be solved
%  language             = programming language of functions file
%  k                    = Search counter
%  field                = cell array of field names of the gui.Term substructure
%*******************************************************************************
function [Problem,Options] = loadMADS(gui,gui_var)

% Transfer Optimization Problem data from GUI into MADS input variables
Problem.nameCache   = gui_var.nameCache;
Problem.typeProblem = gui_var.typeProblem;
nameProblem = get(gui.problem,'String');

ext = fieldnames(gui_var.FileExt);
for k = 1:length(ext)
   Problem.File.(ext{k}) = [nameProblem,gui_var.FileExt.(ext{k})];
end

% Compile non-Matlab Functions file, if possible
Problem.fType = upper(gui_var.problemExt(2));
if ~strcmp(Problem.fType,'M')
   language = nomadm_compile(Problem.fType,Problem.File.F,gui_var.problemExt);
   set(gui.runStatus,'String',['Compiling ',language,' function file']);
end

% Transfer user options from GUI into MADS input variables
Options        = gui_var.Options;
Options.Search = Options.Search(1:Options.nSearches);
[searchTypes]  = getSearchTypesAndLabels(gui,gui_var);
for k = 1:Options.nSearches
   Options.Search(k).type = searchTypes{gui_var.Choice.search(k)};
end
Options.degeneracyScheme = gui_var.Types.degeneracyScheme{[1;2;3;4]'* ...
                           strcmp(get(gui.degeneracyMenu, 'checked'),'on')};
Options.loadCache        = strcmp(get(gui.CacheMenu.load, 'checked'),'on');
Options.countCache       = strcmp(get(gui.CacheMenu.count,'checked'),'on');
Options.countNotInX      = strcmp(get(gui.CacheMenu.countNotInX,'checked'),'on');
Options.removeRedundancy = strcmp(get(gui.OptionsMenu.removeRedundancy,...
                                     'checked'),'on');
Options.runStochastic    = strcmp(get(gui.OptionsMenu.runStochastic, ...
                                     'checked'),'on');
Options.fixCategorical   = strcmp(get(gui.OptionsMenu.fixCategorical, ...
                                     'checked'),'on');
Options.accelerate  = strcmp(get(gui.OptionsMenu.accelerate,  'checked'),'on');
Options.plotFilter  = strcmp(get(gui.OptionsMenu.plotFilter,  'checked'),'on');
Options.saveHistory = strcmp(get(gui.OptionsMenu.saveHistory, 'checked'),'on');
Options.TermFlag.relative = strcmp(get(gui.OptionsMenu.TermFlag.relative, ...
                                      'checked'),'on');
Options.useFilter   = [0;1;2]'*strcmp(get(gui.filterMenu,      'checked'),'on');
Options.useMultiFi  = [0;1;2]'*strcmp(get(gui.multiFiMenu,     'checked'),'on');
Options.scale       = [0;2;10]'*strcmp(get(gui.scaleMenu,      'checked'),'on');
Options.dimSensor   = [0;1;2;3]'*strcmp(get(gui.sensorMenu,    'checked'),'on');
Options.debug       = [0;1;2;3]'*strcmp(get(gui.debugMenu,     'checked'),'on');
Options.plotHistory = [0;1;2]'*strcmp(get(gui.plotHistoryMenu, 'checked'),'on');
Options.runOneIteration  = strcmp(get(gui.RunMenu.oneIteration,'checked'),'on');
Options.runUntilFeasible = strcmp(get(gui.RunMenu.execFeasible,'checked'),'on');
Options.plotColor        = gui_var.Types.plotColors(gui_var.runCount+1);
Options.hplothandle      = gui.axesHistory;
Options.fplothandle      = gui.axesFilter;
Options.stophandle       = gui.stopRun;
Options.delta0           = str2double(get(gui.delta0,      'String'));
Options.deltaMax         = str2double(get(gui.deltaMax,    'String'));
Options.meshRefine       = str2double(get(gui.meshRefine,  'String'));
Options.meshCoarsen      = str2double(get(gui.meshCoarsen, 'String'));
Options.tolCache         = str2double(get(gui.tolCache,    'String'));
Options.hmin             = str2double(get(gui.hmin,        'String'));
Options.hmax             = str2double(get(gui.hmax,        'String'));
Options.hRho             = str2double(get(gui.hRho,        'String'));
Options.EPoll.fTrigger   = str2double(get(gui.ePollXiF,    'String'));
Options.EPoll.hTrigger   = str2double(get(gui.ePollXiH,    'String'));
Options.nRuns            = str2double(get(gui.nRuns,       'String'));

% Process Termination Criteria
field = fieldnames(gui.TermFlag);
for k = 1:length(field)
   Options.Term.(field{k})     = str2double(get(gui.Term.(field{k}),'String'));
   Options.TermFlag.(field{k}) = get(gui.TermFlag.(field{k}),       'Value');
   if ~Options.TermFlag.(field{k}), Options.Term.(field{k}) = Inf; end
end
if Options.TermFlag.relative
   Options.Term.delta = Options.Term.delta*Options.delta0;
end

end

%*******************************************************************************
% edit_gui: Edit parameters.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: MADS-->Edit Parameters)
% Calls:     loadParameters
% VARIABLES:
%  guiEdit  = structure for the Edit Parameters GUI and its object handles
%    .fig   =   handle for the Edit Parameters figure window
%  gui_var  = structure of GUI variables
%  gui.func = handles of all NOMADm callback functions
%  nFields  = number of data fields to be displayed on Edit screen
%  param    = User data storage
%  header   = title of input dialog box
%  labels   = text labels in input dialog box
%  handles  = handles to appropriate GUI figure objects
%  defAns   = default answers that appear in input dialog box
%  newAns   = new answers that appear in input dialog box
%*******************************************************************************
function edit_gui(~,~) %#ok

% Retrieve data on which set of parameters is to be edited
param    = get(gcbo,'UserData');
[labels, handles] = deal(param{:});
defAns   = get(handles, 'String');
if length(defAns) == 1, defAns = {defAns}; end
nFields  = length(labels);
gui      = getappdata(gcbf,'gui');

% Set up parameters used in forming figure window
switch nFields
case {4}
   pos       = 0.585;
   winSize   = 0.28;
   hgt       = [.130, .10];
   rate      = .2115;
   offset    = [1.028, 1.0];
   btnYPos   = .018;
   btnhgt    = .11;
   btnboxhgt = .12;
case {5}
   pos       = 0.575;
   winSize   = 0.340;
   hgt       = [.105, .08];
   rate      = .17;
   offset    = [1.01, .99];
   btnYPos   = .016;
   btnhgt    = .09;
   btnboxhgt = .10;
case {6}
   pos       = 0.575;
   winSize   = 0.55;
   hgt       = [.080, .06];
   rate      = .140;
   offset    = [1.01, .99];
   btnYPos   = .017;
   btnhgt    = .07;
   btnboxhgt = .08;
case {7}
   pos       = 0.575;
   winSize   = 0.55;
   hgt       = [.080, .06];
   rate      = .125;
   offset    = [1.01, .99];
   btnYPos   = .017;
   btnhgt    = .07;
   btnboxhgt = .08;
otherwise
   msg = ['Code for editing ',int2str(nFields),' parameters not available'];
   error('nomadm:edit_gui', msg);
end

%Set up "Edit ..." figure window
guiEdit.fig = figure(...
   'Name',                      get(gcbo,'Label'), ...
   'Tag',                       'Edit_GUI',        ...
   'DefaultUIControlUnits',     'normalized',      ...
   'DefaultUIControlFontUnits', 'normalized',      ...
   'DefaultUIControlFontName',  'Helvetica',       ...
   'DefaultUIControlFontSize',  0.375,             ...
   'DefaultUIControlStyle',     'text',            ...
   'WindowStyle',               'modal',           ...
   'Resize',                    'off',             ...
   'MenuBar',                   'none',            ...
   'NumberTitle',               'off',             ...
   'Units',                     'normalized',      ...
   'Position',                  [0.362, pos-winSize/2, 0.275, winSize]);

% Set up panel for Top displays
uipanel('Parent',guiEdit.fig,'BorderType','beveledin');

for k = 1:nFields
   uicontrol(guiEdit.fig, ...
       'Style',               'text', ...
       'String',              labels{k}, ...
       'HorizontalAlignment', 'left', ...
       'Position',            [.017,offset(1)-k*rate,.973,hgt(1)]);
   guiEdit.newAns(k) = uicontrol(guiEdit.fig, ...
       'Style',               'edit', ...
       'String',              defAns{k}, ...
       'BackgroundColor',     'white', ...
       'HorizontalAlignment', 'left', ...
       'FontSize',            .50, ...
       'Position',            [.017,offset(2)-k*rate,.973,hgt(2)]);
end

% The Done and Cancel Buttons
uicontrol(guiEdit.fig, 'Style','frame', 'ForegroundColor','black', ...
   'Position',[.575, .0128, .19, btnboxhgt]);
guiEdit.done = uicontrol(guiEdit.fig,         ...
   'Style',      'pushbutton',                ...
   'String',     'OK',                        ...
   'Fontsize',   .45,                         ...
   'Position',   [.58, btnYPos, .18, btnhgt], ...
   'UserData',   {handles},                   ...
   'Callback',   gui.func.loadEdit);
guiEdit.cancel = uicontrol(guiEdit.fig,       ...
   'Style',      'pushbutton',                ...
   'String',     'Cancel',                    ...
   'Fontsize',   .45,                         ...
   'Position',   [.78, btnYPos, .19, btnhgt], ...
   'Callback',   'close(gcf)');

setappdata(guiEdit.fig,'guiEdit',guiEdit);
end

%*******************************************************************************
% loadEdit:  Updates screen for edited parameter choices.
% ------------------------------------------------------------------------------
% Called by: edit_gui
% VARIABLES:
%  guiEdit   = structure for the Edit Parameters GUI and its object handles
%    .newAns = handles for Edit Parameters figure window objects
%  handles   = handles to NOMADm GUI figure objects
%*******************************************************************************
function loadEdit(~,~) %#ok

guiEdit  = getappdata(findobj('Tag','Edit_GUI'),'guiEdit');
handles  = get(gcbo,'UserData');
set(handles{:}, {'String'}, get(guiEdit.newAns,'String'));
close(gcf);
end

%*******************************************************************************
% search_gui:  Displays a user input screen to set Search parameters.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: MADS-->Select Search Strategies)
% Calls:     modifySearchScreen, loadUserSearchFile, loadSearchOptions
% VARIABLES:
%  guiSearch         = structure for the Search GUI and its object handles
%    .fig            =   handle for Search figure window
%    .n              =   handle for Number of Search Types popup menu
%    .SurOptimizer   =   handle for surrogate optimizer popup menu
%    .Screen         = structure of handles for each k of n Searches
%      .label        =   handle for Search text label
%      .type         =   handle for popupmenu of possible Search types
%      .nIter        =   handle for number of iterations field
%      .nPoints      =   handle for number of Search points field
%      .file         =   handle for string name of user Search file
%      .sfile        =   handle for string names of surrogate files
%    .done           =   handle for Done pushbutton
%    .cancel         =   handle for Cancel pushbutton   
%  gui_var           = structure of GUI variables
%    .maxSearches    = maximum number of Search types that can be selected
%    .Labels         = long text labels used in popup menus
%      .Search       =   labels used in each Search type popup menu
%      .optimizer    =   labels for the surrogate optimizer to be used
%    .Choice.Search  = integers recording Search type popup menu choices
%    .Options.Search = vector of structures for each k of n Searches
%      .nIter        =   number of iterations
%      .nPoints      =   number of Search points
%      .file         =   string containing name of optional Search file
%  gui.func          = handles of all NOMADm callback functions
%  k                 = Search type counter
%  row               = row location of current Search figure window object
%  pathStr           = temporary storage for Search file path
%  filename          = temporary storage for Search filename
%  fileExt           = temporary storage for Search filename extension
%*******************************************************************************
function search_gui(~,~) %#ok

gui_var  = getappdata(gcbf,'gui_var');
gui      = getappdata(gcbf,'gui');

%Set up "Search Options" figure window
guiSearch.fig = figure(...
   'Name',                      'Set Options for MADS SEARCH step', ...
   'Tag',                       'Search_GUI', ...
   'DefaultUIControlUnits',     'normalized', ...
   'DefaultUIControlFontUnits', 'normalized', ...
   'DefaultUIControlFontName',  'Helvetica',  ...
   'DefaultUIControlFontSize',  0.375, ...,
   'DefaultUIControlStyle',     'text', ...
   'WindowStyle',               'modal', ...
   'Units',                     'normalized', ...
   'Position',                  [0.15 0.1 0.8 0.7], ...
   'MenuBar',                   'none', ...
   'NumberTitle',               'off');

% Set up panel for Top displays
panelColor      = get(guiSearch.fig,'DefaultUIControlBackgroundColor');
guiSearch.panel = uipanel('Parent',guiSearch.fig, ...
                          'Position',[0, .9, 1, .1], ...
                          'BorderType','beveledin', ...
                          'BackgroundColor',panelColor);

% Display Number of Search Types and Surrogate Optimizer
uicontrol(guiSearch.panel, ...
    'Style',           'text', ...
    'String',          'Number of Search Types: ', ...
    'Position',        [0 0 .18 .6]);
guiSearch.n = uicontrol(guiSearch.panel, ...
    'Style',           'popupmenu', ...
    'String',          {'0','1','2','3','4','5','6','7','8'}, ...
    'BackgroundColor', 'white', ...
    'Value',           gui_var.Options.nSearches+1, ...
    'Position',        [.18 .03 .05 .65], ...
    'Callback',        gui.func.modifySearchScreen);
uicontrol(guiSearch.panel, ...
    'Style',           'text', ...
    'String',          'Surrogate Optimizer: ', ...
    'Position',        [.24 0 .15 .6]);
guiSearch.SurOptimizer = uicontrol(guiSearch.panel, ...
    'Style',           'popupmenu', ...
    'String',          gui_var.Labels.optimizer, ...
    'BackgroundColor', 'white', ...
    'Value',           gui_var.Choice.optimizer, ...
    'Position',        [.39 .03 .33 .65]);
guiSearch.mvp1Surrogate = uicontrol(guiSearch.panel, ...
    'Style',           'checkbox', ...
    'String',          'Use single surrogate for MVPs', ...
    'Value',           gui_var.Options.mvp1Surrogate, ...
    'Position',        [.75 .03 .23 .65]);

% Text headers for the Search parameters
headerPosX   = {0.05,0.29,0.40,0.49,0.56,0.70,0.85};
headerLength = {0.23,0.10,0.08,0.06,0.13,0.14,0.14};
headerText   = {'SearchStrategy','# Iterations', '# Points','Complete', ...
                'Surrogate Files','User File','Apply Surrogate to:'};
headerFld    = {'','','','complete','sfile','file','param'};
for k = 1:3
   uicontrol(guiSearch.fig,'Style','text','String',headerText{k}, ...
             'Position', [headerPosX{k}, .80, headerLength{k}, .06]);
end
for k = 4:7
   guiSearch.Labels.(headerFld{k}) = uicontrol(guiSearch.fig,'Style','text', ...
             'String',    headerText{k}, ...
             'Position', [headerPosX{k}, .80, headerLength{k}, .06]);
end

% Disable Nadaraya-Watson option unless the problem is stochastic
[~,searchLabels] = getSearchTypesAndLabels(gui,gui_var);

% Main loop for each of n Searches
for k = 1:gui_var.maxSearches
   row = .76 - .06*(k-1);

   guiSearch.Screen(k).label = uicontrol(guiSearch.fig, ...
      'Style',           'text', ...
      'String',          [int2str(k), ':'], ...
      'Position',        [.01, row-.0075, .03, .06]);

   % The data fields for each Search
   guiSearch.Screen(k).type = uicontrol(guiSearch.fig, ...
      'Style',           'popupmenu', ...
      'String',          searchLabels, ...
      'BackgroundColor', 'white', ...
      'Value',           gui_var.Choice.search(k), ...
      'Position',        [headerPosX{1}, row, headerLength{1}, .06], ...
      'UserData',        {k,gui_var.Choice.search(k), ...
                         gui_var.Options.Search(k).local, ...
                         gui_var.Options.Search(k).merit, ...
                         gui_var.Options.Search(k).cblgs.nGoal}, ...
      'Callback',        gui.func.loadUserSearchFile);

   guiSearch.Screen(k).nIter = uicontrol(guiSearch.fig, ...
      'Style',           'edit', ...
      'String',          int2str(gui_var.Options.Search(k).nIter), ...
      'BackgroundColor', 'white', ...
      'FontSize',        .64, ...
      'Position',        [headerPosX{2}, row+.02, headerLength{2}, .04]);
   guiSearch.Screen(k).nPoints = uicontrol(guiSearch.fig, ...
      'Style',           'edit', ...
      'String',          int2str(gui_var.Options.Search(k).nPoints), ...
      'BackgroundColor', 'white', ...
      'FontSize',        .64, ...
      'Position',        [headerPosX{3}, row+.02, headerLength{3}, .04]);
   guiSearch.Screen(k).complete = uicontrol(guiSearch.fig, ...
      'Style',           'checkbox', ...
      'Value',           gui_var.Options.Search(k).complete, ...
      'Position',        [ headerPosX{4}+0.025, row+.02, .03, .04]);
   
   [pathStr,fileName1,fileExt1] = fileparts(gui_var.Options.Search(k).sfile{1});
   [~,      fileName2,fileExt2] = fileparts(gui_var.Options.Search(k).sfile{2});
   searchLabel = [fileName1,fileExt1,'/',fileName2,fileExt2];
   if strcmp(searchLabel,'/'), searchLabel = ''; end
   guiSearch.Screen(k).sfile  = uicontrol(guiSearch.fig, ...
      'Style',           'text', ...
      'String',          searchLabel, ...
      'ForegroundColor', 'red', ...
      'Position',        [headerPosX{5}, row, headerLength{5}, .06], ...
      'UserData',        pathStr);
   if isempty(gui_var.Options.Search(k).file)
      sFile = '';
   else
      sFile = func2str(gui_var.Options.Search(k).file);
   end
%    [pathStr,filename,fileExt] = fileparts(gui_var.Options.Search(k).file);
   guiSearch.Screen(k).file   = uicontrol(guiSearch.fig, ...
      'Style',           'text', ...
      'String',          sFile, ...
      'ForegroundColor', 'red', ...
      'Position',        [headerPosX{6}, row, headerLength{6}, .06], ...
      'UserData',        gui_var.Options.Search(k).file);
   guiSearch.Screen(k).param   = uicontrol(guiSearch.fig, ...
      'Style',           'popupmenu', ...
      'String',          {'Objective Function', 'Other Data'}, ...
      'Value',           gui_var.Options.Search(k).param + 1, ...
      'BackgroundColor', 'white', ...
      'Position',        [headerPosX{7}, row, headerLength{7}, .06]);
end

% The Done and Cancel Buttons
guiSearch.done = uicontrol(guiSearch.fig, ...
   'Style',      'pushbutton', ...
   'String',     'Done', ...
   'FontWeight', 'bold', ...
   'Position',   [.33, .04, .10, .08], ...
   'Callback',   gui.func.loadSearchOptions);
guiSearch.cancel = uicontrol(guiSearch.fig, ...
   'Style',      'pushbutton', ...
   'String',     'Cancel', ...
   'FontWeight', 'bold', ...
   'Position',   [.50, .04, .10, .08], ...
   'Callback',   'close(gcf)');

setappdata(guiSearch.fig,'guiSearch',guiSearch);
modifySearchScreen(2,0);
end

%*******************************************************************************
% getSearchTypesAndLabels: Delete stochastic options if search is not stochastic
% ------------------------------------------------------------------------------
% Called by: search_gui
% VARIABLES:
%  searchTypes        = cell array of Search types without stochastic options
%  searchLabels       = cell array of Search labels without stochastic options
%  gui                = structure of all GUI object handles
%    .OptionsMenu     =   structure of Options menu objects
%      .runStochastic =     handle to flag for running as stochastic problem
%  gui_var            = structure of all GUI variables
%    .Labels.search   =   cell array of text labels for all Search options
%    .Types.search    =   cell array of Search type strings
%  ind                = indices of search types and labels that are stochastic
%*******************************************************************************
function [searchTypes,searchLabels] = getSearchTypesAndLabels(gui,gui_var)

searchTypes  = gui_var.Types.search;
searchLabels = gui_var.Labels.search;
if ~strcmp(get(gui.OptionsMenu.runStochastic,'checked'),'on')
   ind = find(contains(searchTypes,'NW'));
   searchTypes(ind)  = [];
   searchLabels(ind) = [];
end
end

%*******************************************************************************
% modifySearchScreen: Modify the Search Input Screen.
% ------------------------------------------------------------------------------
% Called by: search_gui (callback: Number of Search Types popupmenu)
%            search_gui (directly), loadSurrogateOptions, loadUserSearchFile
% VARIABLES:
%  fig            = temporary storage of figure handles
%  gui_var        = structure of all GUI variables
%    .nSearches   =   handle for Number of Search Types popup menu
%    .maxSearches =   maximum number of Search types allowed
%  guiSearch      = structure of all Search window object handles
%    .n           =   number of Search types popup menu
%    .Screen(k)   =   Search fields for the Search #k
%      .label     =     text label appearing on the Search Screen
%      .type      =     type of Search
%      .nIter     =     number of iterations
%      .nPoints   =     number of Search points
%      .file      =     file used during Search
%    .done        =   handle for Search input screen Done button
%    .cancel      =   handle for Search Input Screen Cancel button   
%  nSearches      = number of Search types selected
%  isSurrogate    = flags to indicate if Searches use surrogates
%  handles        = handles of objects whose visibility is toggled
%  pos            = screen position of the final Search field
%*******************************************************************************
function modifySearchScreen(~,~)

% Get application data
onoff     = {'on','off'};
fig       = findobj('Tag','NOMADm_GUI');
gui_var   = getappdata(fig,'gui_var');
fig       = findobj('Tag','Search_GUI');
guiSearch = getappdata(fig,'guiSearch');

% Turn on/off appropriate Search fields
nSearches   = get(guiSearch.n,'Value') - 1;
isSurrogate = zeros(1,nSearches);
hasUFile    = zeros(1,nSearches);
hasSFile    = zeros(1,nSearches);
for k = 1:gui_var.maxSearches
   if k <= nSearches
      isSurrogate(k) = ~isempty(strfind(gui_var.Labels.search{ ...
                        get(guiSearch.Screen(k).type,'Value')},'Surrogate'));
      hasUFile(k)   = ~isempty(deblank(get(guiSearch.Screen(k).file,'String')));
      searchFileStr = regexprep(get(guiSearch.Screen(k).sfile,'String'),'/','');
      hasSFile(k) = ~isempty(deblank(searchFileStr));
      set(guiSearch.Screen(k).label,    'Visible','on');
      set(guiSearch.Screen(k).type,     'Visible','on');
      set(guiSearch.Screen(k).nIter,    'Visible','on');
      set(guiSearch.Screen(k).nPoints,  'Visible','on');
      set(guiSearch.Screen(k).complete, 'Visible','on');
      set(guiSearch.Screen(k).sfile,    'Visible','on');
      set(guiSearch.Screen(k).file,     'Visible','on');
   else
      set(guiSearch.Screen(k).label,    'Visible','off');
      set(guiSearch.Screen(k).type,     'Visible','off','Value', 1);
      set(guiSearch.Screen(k).nIter,    'Visible','off','String','1');
      set(guiSearch.Screen(k).nPoints,  'Visible','off','String','1');
      set(guiSearch.Screen(k).complete, 'Visible','off');
      set(guiSearch.Screen(k).sfile,    'Visible','off','String','');
      set(guiSearch.Screen(k).file,     'Visible','off','String','');
      set(guiSearch.Screen(k).param,    'Visible','off','Value',1);
   end
end

% Make appropriate fields visible, as appropriate
shandles = [guiSearch.Labels.sfile; [guiSearch.Screen(1:nSearches).sfile]'];
uhandles = [guiSearch.Labels.file;  [guiSearch.Screen(1:nSearches).file]'];

set(shandles,'Visible',onoff{2-any(hasSFile(1:nSearches))});
set(uhandles,'Visible',onoff{2-any(hasUFile(1:nSearches))});
set(guiSearch.Labels.param,'Visible',onoff{2-any(isSurrogate(1:nSearches))});
for k = 1:nSearches
   set([guiSearch.Screen(k).param]','Visible',onoff{2-isSurrogate(k)});
end

% Move Done and Cancel buttons
if nSearches
   pos = get(guiSearch.Screen(nSearches).type,'Position');
else
   pos = [.01, .82];
end
set(guiSearch.done,   'Position', [.33, pos(2)-.16, .10, .08]);
set(guiSearch.cancel, 'Position', [.50, pos(2)-.16, .10, .08]);
end

%*******************************************************************************
% loadUserSearchFile: Load a user-specified Search file.
% ------------------------------------------------------------------------------
% Called by: search_gui (callback: each Search Strategy popup menu)
% Calls:     cmaes_gui, dace_gui, nw_gui, rbf_gui, modifySearchScreen
% VARIABLES:
%  fig                   = temporary storage of figure handles
%  gui_var               = structure of all GUI variables
%    .Types.search       =   list of possible Search types
%    .maxSearches        =   maximum number of Searches that can be done
%  guiSearch             = handles for all Search window objects
%    .Screen.file        =   handle for Search figure window user file
%  newValue              = number associated with Search menu choice
%  userData              = temporary storage of GUI object user data
%  k                     = Search number
%  previousType          = previously selected Search type
%  gui.func              = handles of all NOMADm callback functions
%    .loadUserSearchFile = handle for the function of the same name
%  SName, SPath          = name and path of user Search file
%  spec                  = file types for use in input dialog boxes
%*******************************************************************************
function loadUserSearchFile(~,~) %#ok

spec = {'*.m',        'Matlab M-files (*.m)'; ...
        '*.f; *.for', 'Fortran files (*.f,*.for)'; ...
        '*.c; *.C',   'C/C++ files (*.c,*.C)'};

% Get application data
fig         = findobj('Tag','NOMADm_GUI');
gui_var     = getappdata(fig,'gui_var');
gui         = getappdata(fig,'gui');
fig         = findobj('Tag','Search_GUI');
guiSearch   = getappdata(fig,'guiSearch');
newValue    = get(gcbo, 'Value');
userData    = get(gcbo, 'UserData');
[searchTypes] = getSearchTypesAndLabels(gui,gui_var);
[k,previousType,local,merit,nGoal] = deal(userData{:});

% Choose action based on which search type was chosen
switch searchTypes{newValue}

   case {'None'}
      set(guiSearch.Screen(k).nPoints,  'String','1');
      set(guiSearch.Screen(k).complete, 'Value',0);
      set(guiSearch.Screen(k).file,     'String','');
      set(guiSearch.Screen(k).sfile,    'String','');    
      set(gcbo,'UserData',{k,newValue,local,merit,nGoal}); 

      nSearches = get(guiSearch.n,'Value') - 1;     
      if k == nSearches
         set(guiSearch.n,'Value',get(guiSearch.n,'Value')-1);
      end

   % Load a custom Search or Surrogate file 
   case {'Custom','CustomS'}
      sName = uigetfile(spec,['Choose File for Search #',int2str(k)]);
      if (sName)
         sFunc = str2func(strtok(sName,'.'));
         set(guiSearch.Screen(k).sfile,'String','');
         set(guiSearch.Screen(k).file, 'UserData', sFunc);
         set(guiSearch.Screen(k).file, 'String',   sName);
         set(gcbo, 'UserData',{k,newValue,local,merit,nGoal});
      else
         set(gcbo, 'Value', previousType);
      end

   % Load the CMA-ES GUI window
   case {'CMAES'}
      cmaes_gui(k,gui_var,gui);

   % Load a DACE surrogate
   case {'DACE'}
      dace_gui(k,gui_var,gui);

   % Load a NW surrogate
   case {'NW'}
      nw_gui(k,gui_var,gui);

   % Load a NW surrogate
   case {'RBF'}
      rbf_gui(k,gui_var,gui);

   % Load an SPS-DACE surrogate
   case {'SPS-DACE'}
      sName = uigetfile(spec,['Choose File for Search #',int2str(k)]);
      if sName
         guiSurrogate = dace_gui(k,gui_var,gui);
         sFunc = str2func(strtok(sName,'.'));
         set(guiSurrogate.fig,'UserData',{sFunc,sName});
      else
         set(gcbo,'Value',previousType);
      end

   % Load an SPS-NW surrogate
   case {'SPS-NW'}
      sName = uigetfile(spec,['Choose File for Search #',int2str(k)]);
      if sName
         guiSurrogate = nw_gui(k,gui_var,gui);
         sFunc = str2func(strtok(sName,'.'));
         set(guiSurrogate.fig,'UserData',{sFunc,sName});
      else
         set(gcbo,'Value',previousType);
      end

   % Load an SPS-NW surrogate
   case {'SPS-RBF'}
      sName = uigetfile(spec,['Choose File for Search #',int2str(k)]);
      if sName
         guiSurrogate = rbf_gui(k,gui_var,gui);
         sFunc = str2func(strtok(sName,'.'));
         set(guiSurrogate.fig,'UserData',{sFunc,sName});
      else
         set(gcbo,'Value',previousType);
      end

   % Do not load a file
   otherwise
      set(guiSearch.Screen(k).file, 'String','');
      set(guiSearch.Screen(k).sfile,'String','');
      set(gcbo,'UserData',{k,newValue,local,merit,nGoal});
end

modifySearchScreen(2,0);
end

%*******************************************************************************
% loadSearchOptions: Load selected Search parameters.
% ------------------------------------------------------------------------------
% Called by: search_gui (callback: Done and Cancel pushbuttons)
% Calls:     updateSearchLabels
% VARIABLES:
%  k                = Search number
%  fig              = temporary storage of the appropriate figure handle
%  gui_var          = structure of all GUI handles and variables
%    .Choice        =   structure of user choices
%      .search(k)   =     user Search choices
%      .optimizer   =     choice of surrogate optimizer
%    .Options       =   structure of MADS parameters settings
%      .Search(k)   =     structure of user-selected k-th Search
%      .nSearches   =     number of Search types to be used
%    .Types         =   lists of possible types
%      .Search      =     list of possible Search types
%      .optimizer   =     list of possible surrogate optimizers
%    .noGrad        =   flag indicating initial availability of gradients
%  guiSearch        = handles for all Search window objects
%    .n             =   handle for the Number of Searches field
%    .SurOptimizer  =   handle for the Surrogate Optimizer popup menu
%    .Screen        =   object handles for each Search
%      .type        =     string identifying Search type
%      .label       =     long text label for Search type
%      .nIter       =     number of iterations to perform Search
%      .nPoints     =     number of Search points
%      .file        =     optional user file defining Search
%  loadSearch       = flag indicating if Search options will be loaded
%  Search           = temporary storage of gui_var.Options.Search(k)
%  maxDisplay       = number of Search Types displayed on the main GUI
%  searchLabel      = Search label that appears on the main GUI
%*******************************************************************************
function loadSearchOptions(~,~) %#ok

fig       = findobj('Tag','NOMADm_GUI');
gui       = getappdata(fig,'gui');
gui_var   = getappdata(fig,'gui_var');
guiSearch = getappdata(gcbf,'guiSearch');

[searchTypes,searchLabels]   = getSearchTypesAndLabels(gui,gui_var);  
gui_var.Options.nSearches    = get(guiSearch.n,'Value') - 1;
gui_var.Choice.optimizer     = get(guiSearch.SurOptimizer,'Value');
gui_var.Options.SurOptimizer = gui_var.Types.optimizer{gui_var.Choice.optimizer};
gui_var.Options.mvp1Surrogate = get(guiSearch.mvp1Surrogate,'Value');
for k = 1:gui_var.Options.nSearches
   gui_var.Choice.search(k) = get(guiSearch.Screen(k).type, 'Value');
   Search.type  = searchTypes{gui_var.Choice.search(k)};
   Search.label = searchLabels{gui_var.Choice.search(k)};
   if (gui_var.noGrad && strcmp(Search.type, 'GPollI'))
      uiwait(msgbox('No derivatives available for this problem', ...
                   ['Error in Search Type #', int2str(k)],'error','modal'));
      return
   end
   Search.nIter    = str2double(get(guiSearch.Screen(k).nIter,  'String'));
   Search.nPoints  = str2double(get(guiSearch.Screen(k).nPoints,'String'));
   Search.complete = get(guiSearch.Screen(k).complete,'Value');
   pathStr         = get(guiSearch.Screen(k).sfile,'UserData');
   searchLabel     = get(guiSearch.Screen(k).sfile,'String');
   if contains(searchLabel,'/')
      sfile = split(searchLabel,'/');
   else
      sfile = {searchLabel,''};
   end
   Search.sfile       = fullfile(pathStr,sfile);
   Search.file        = get(guiSearch.Screen(k).file, 'UserData');
   UserData           = get(guiSearch.Screen(k).type, 'UserData');
   Search.local       = UserData{3};
   Search.merit       = UserData{4};
   Search.param       = get(guiSearch.Screen(k).param, 'Value') - 1;
   Search.cblgs.nGoal = UserData{5};
   gui_var.Options.Search(k) = Search;
end

% Update GUI Search fields
maxDisplay  = length(gui.searchLabel);
searchLabel = updateSearchLabels(maxDisplay, ...
                    gui_var.Options.nSearches, gui_var.Options.Search);
for k = 1:maxDisplay
   set(gui.searchLabel(k), 'String', searchLabel{k});
end

setappdata(fig,'gui_var',gui_var);
close(guiSearch.fig);
end

%*******************************************************************************
% poll_gui:  Displays a user input screen to set Poll parameters.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: MADS-->Select Poll Parameters)
% Calls:     loadPollOptions
% VARIABLES:
%  guiPoll           = structure for the Search GUI and its object handles
%    .fig            =   handle for Search figure window
%    .Poll           =   handles for Poll parameters
%      .label        =     Poll text label
%      .type         =     Popupmenu of possible Poll types
%      .order        =     Popupmenu of possible Poll ordering types
%      .center       =     Index number of Poll center
%      .complete     =     checkbox indicating complete poll
%    .EPoll          =   handles for Extended Poll parameters
%      .completeN    =     checkbox indicating complete discrete neighbor poll
%      .label        =     Poll text label
%      .type         =     Popupmenu of possible Poll types
%      .order        =     Popupmenu of possible Poll ordering types
%      .center       =     Index number of Poll center
%      .complete     =     checkbox indicating complete poll
%    .done           =   handle for Done pushbutton
%    .cancel         =   handle for Cancel pushbutton   
%  gui_var           = structure of GUI variables
%    .Labels         = long text labels used in popup menus
%      .Poll         =   labels used in each Search type popup menu
%    .Choice.Poll    = integers recording Poll type popup menu choices
%    .Options.Poll   = structures for Poll step data
%  gui.func          = handles of all NOMADm callback functions
%*******************************************************************************
function poll_gui(~,~) %#ok

gui_var  = getappdata(gcbf,'gui_var');
gui      = getappdata(gcbf,'gui');

% Limit choices if no derivatives are available
PollLabels  = gui_var.Labels.poll;
if gui_var.noGrad
   PollLabels(strncmp(PollLabels,'Gradient',8)) = [];
end

% Reset Poll choice index if it exceeds the modified Poll type list length
PollChoice = gui_var.Choice.Poll.type;
if PollChoice > length(PollLabels)
   PollChoice = 1;
end
EPollChoice = gui_var.Choice.EPoll.type;
if EPollChoice > length(PollLabels)
   EPollChoice = 1;
end

%Set up "Poll Options" figure window
guiPoll.fig = figure(...
   'Name',                      'Set Options for MADS POLL step', ...
   'Tag',                       'Poll_GUI', ...
   'DefaultUIControlUnits',     'normalized', ...
   'DefaultUIControlFontUnits', 'normalized', ...
   'DefaultUIControlFontName',  'Helvetica',  ...
   'DefaultUIControlFontSize',  0.4, ...,
   'DefaultUIControlStyle',     'text', ...
   'WindowStyle',               'modal', ...
   'Units',                     'normalized', ...
   'Position',                  [0.12 0.15 0.8 0.7], ...
   'MenuBar',                   'none', ...
   'NumberTitle',               'off');

%Set up panel for Top displays
uipanel('BorderType','beveledin');

% Position references
tab1 = .02;
tab2 = .18;
tab3 = .59;
thgt = .09;
plen = .38;

% Text headers for the POLL parameters
uicontrol(guiPoll.fig, 'Style', 'text', ...
    'String',   'POLL Step', ...
    'Position', [tab2, .85, plen, thgt]);
uicontrol(guiPoll.fig, 'Style', 'text', ...
    'String',   'EXTENDED POLL Step', ...
    'Position', [tab3, .85, plen, thgt]);

uicontrol(guiPoll.fig, 'Style', 'text', ...
    'String',              'Poll Directions:', ...
    'HorizontalAlignment', 'left', ...
    'Position',            [tab1, .75, .15, thgt]);
guiPoll.type = uicontrol(guiPoll.fig, ...
   'Style',           'popupmenu', ...
   'String',          PollLabels, ...
   'BackgroundColor', 'white', ...
   'Value',           PollChoice, ...
   'Position',        [tab2, .75, plen, thgt]);
guiPoll.Etype = uicontrol(guiPoll.fig, ...
   'Style',           'popupmenu', ...
   'String',          PollLabels, ...
   'BackgroundColor', 'white', ...
   'Value',           EPollChoice, ...
   'Position',        [tab3, .75, plen, thgt]);

uicontrol(guiPoll.fig, 'Style', 'text', ...
    'String',              'Poll Order:', ...
    'HorizontalAlignment', 'left', ...
    'Position',            [tab1, .65, .15, thgt]);
guiPoll.order = uicontrol(guiPoll.fig, ...
   'Style',           'popupmenu', ...
   'String',          gui_var.Labels.pollOrder, ...
   'BackgroundColor', 'white', ...
   'Value',           gui_var.Choice.Poll.order, ...
   'Position',        [tab2, .65, plen, thgt]);
guiPoll.Eorder = uicontrol(guiPoll.fig, ...
   'Style',           'popupmenu', ...
   'String',          gui_var.Labels.pollOrder, ...
   'BackgroundColor', 'white', ...
   'Value',           gui_var.Choice.EPoll.order, ...
   'Position',        [tab3, .65, plen, thgt]);

uicontrol(guiPoll.fig, 'Style', 'text', ...
    'String',              'Poll Center:', ...
    'HorizontalAlignment', 'left', ...
    'Position',            [tab1, .55, .15, thgt]);
guiPoll.center = uicontrol(guiPoll.fig, ...
   'Style',           'popupmenu', ...
   'String',          gui_var.Labels.pollCenter, ...
   'BackgroundColor', 'white', ...
   'Value',           gui_var.Choice.Poll.center, ...
   'Position',        [tab2, .55, plen, thgt]);
guiPoll.Ecenter = uicontrol(guiPoll.fig, ...
   'Style',           'popupmenu', ...
   'String',          gui_var.Labels.pollCenter, ...
   'BackgroundColor', 'white', ...
   'Value',           gui_var.Choice.EPoll.center, ...
   'Position',        [tab3, .55, plen, thgt]);

uicontrol(guiPoll.fig, 'Style', 'text', ...
    'String',              'Complete Polling', ...
    'HorizontalAlignment', 'left', ...
    'Position',            [.40, .40, .30, .10]);
guiPoll.complete = uicontrol(guiPoll.fig, ...
   'Style',           'checkbox', ...
   'String',          ' Poll step', ...
   'Value',           gui_var.Options.Poll.complete, ...
   'Position',        [.40, .35, .30, thgt]);
guiPoll.Ncomplete = uicontrol(guiPoll.fig, ...
   'Style',           'checkbox', ...
   'String',          ' Discrete Neighbor Poll', ...
   'Value',           gui_var.Options.EPoll.completeN, ...
   'Position',        [.40, .28, .30, thgt]);
guiPoll.Ecomplete = uicontrol(guiPoll.fig, ...
   'Style',           'checkbox', ...
   'String',          ' Extended Poll step', ...
   'Value',           gui_var.Options.EPoll.complete, ...
   'Position',        [.40, .21, .30, thgt]);

% The Done and Cancel Buttons
guiPoll.done = uicontrol(guiPoll.fig, ...
   'Style',      'pushbutton', ...
   'String',     'Done', ...
   'FontWeight', 'bold', ...
   'Position',   [.36, .04, .10, thgt], ...
   'Callback',   gui.func.loadPollOptions);
guiPoll.cancel = uicontrol(guiPoll.fig, ...
   'Style',      'pushbutton', ...
   'String',     'Cancel', ...
   'FontWeight', 'bold', ...
   'Position',   [.57, .04, .10, thgt], ...
   'Callback',   'close(gcf)');

setappdata(guiPoll.fig,'guiPoll',guiPoll);
end

%*******************************************************************************
% loadPollOptions: Load selected Poll parameters.
% ------------------------------------------------------------------------------
% Called by: poll_gui (callback: Done and Cancel pushbuttons)
% Calls:     updatePollLabels
% VARIABLES:
%  fig              = temporary storage of the appropriate figure handle
%  gui_var          = structure of all GUI handles and variables
%    .Choice        =   structure of user choices
%      .Poll        =     user Poll choices
%    .Options       =   structure of MADS parameters settings
%      .Poll        =     structure of user-selected Poll
%    .Types         =   lists of possible types
%      .Poll        =     list of possible Poll types
%  guiPoll           = structure for the Poll GUI and its object handles
%    .fig            =   handle for Poll figure window
%    .Poll           =   handles for Poll parameters
%      .label        =     Poll text label
%      .type         =     Popupmenu of possible Poll types
%      .order        =     Popupmenu of possible Poll ordering types
%      .center       =     Index number of Poll center
%      .complete     =     checkbox indicating complete poll
%    .EPoll          =   handles for Extended Poll parameters
%      .completeN    =     checkbox indicating complete discrete neighbor poll
%      .label        =     Poll text label
%      .type         =     Popupmenu of possible Poll types
%      .order        =     Popupmenu of possible Poll ordering types
%      .center       =     Index number of Poll center
%      .complete     =     checkbox indicating complete poll
%    .done           =   handle for Done pushbutton
%    .cancel         =   handle for Cancel pushbutton   
%*******************************************************************************
function loadPollOptions(~,~) %#ok

fig     = findobj('Tag','NOMADm_GUI');
gui     = getappdata(fig,'gui');
gui_var = getappdata(fig,'gui_var');
guiPoll = getappdata(gcbf,'guiPoll');

gui_var.Choice.Poll.type    = get(guiPoll.type,    'Value');
gui_var.Choice.Poll.order   = get(guiPoll.order,   'Value');
gui_var.Choice.Poll.center  = get(guiPoll.center,  'Value');
gui_var.Choice.EPoll.type   = get(guiPoll.Etype,   'Value');
gui_var.Choice.EPoll.order  = get(guiPoll.Eorder,  'Value');
gui_var.Choice.EPoll.center = get(guiPoll.Ecenter, 'Value');
gui_var.Options.Poll.type  = gui_var.Types.poll{gui_var.Choice.Poll.type};
gui_var.Options.Poll.order = gui_var.Types.pollOrder{gui_var.Choice.Poll.order};
gui_var.Options.Poll.center = gui_var.Choice.Poll.center - 1;
gui_var.Options.EPoll.type  = gui_var.Types.poll{gui_var.Choice.EPoll.type};
gui_var.Options.EPoll.order = gui_var.Types.pollOrder{gui_var.Choice.EPoll.order};
gui_var.Options.EPoll.center = gui_var.Choice.EPoll.center - 1;
gui_var.Options.Poll.complete   = get(guiPoll.complete,  'Value');
gui_var.Options.EPoll.completeN = get(guiPoll.Ncomplete, 'Value');
gui_var.Options.EPoll.complete  = get(guiPoll.Ecomplete, 'Value');
setappdata(fig,'gui_var',gui_var);

% Update GUI Poll fields
set(gui.pollStrategy,'String',gui_var.Labels.poll{gui_var.Choice.Poll.type});
set(gui.pollOrder,'String',gui_var.Labels.pollOrder{gui_var.Choice.Poll.order});
set(gui.pollCenter,'String',gui_var.Labels.pollCenter{gui_var.Choice.Poll.center});
close(guiPoll.fig);
end

%*******************************************************************************
% cmaes_gui:  Displays a user input screen to set NW Toolbox parameters.
% ------------------------------------------------------------------------------
% Called by: loadUserSearchFile
% Calls:     loadCMAESOptions
% VARIABLES:
%  k                  = Search number for this RBF screen
%  guiCMAES           = structure of all CMAES window object handles
%    .fig             =   CMAES figure window
%    .sigma           =   handle for CMAES sigma field
%    .Stop            =   handles for CMAES stopping criteria 
%      .fitness       =     handle for minimal fitness function value
%      .maxFunEvals   =     handle for maximum number of function evaluations
%      .maxIter       =     handle for maximum number of iterations
%      .tolX          =     handle for minimal variable change
%      .tolUpX        =     handle for maximal variable change
%      .tolFun        =     handle for minimal objective function value change
%    .evalInitialX    =   handle for checkbox to evaluate initial point
%    .incPopSize      =   handle for population size increment factor
%    .popSize         =   handle for initial population size
%    .nParents        =   handle for CMAES number of parents
%    .nRestarts       =   handle for number of restarts
%    .recombWeights   =   handle for recombination weighting strategy popupmenu
%    .done            =   the Done pushbutton
%    .cancel          =   the Cancel pushbutton
%  fig                = handle for the NOMADm figure window
%  gui_var            = structure of all GUI variables
%    .Labels          =   labels for CMAES function popup menus
%      .recombWeights =     labels for CMAES recombination weights
%  gui.func           = structure of all GUI function handles
%  Options            = structure of MADS options
%    .cmaes           =   user-chosen CMAES parameters
%      .sigma         =     standard deviations for each variable
%      .Stop          =     structure of stopping conditions
%        .Fitness     =       fitness value lower threshold
%        .MaxFunEvals =       maximum number of function evaluations
%        .MaxIter     =       maximum number of iterations
%        .TolX        =       minimum distance between x-values
%        .TolUpX      =       maximum distance between x-values
%        .TolFun      =       minimum difference between function values
%      .evalInitialX  =     flag for evaluating the initial point
%      .nRestarts     =     number of restarts
%      .incPopSize    =     factor for multiplying population size
%      .popSize       =     population size 
%      .nParents      =     number of parents
%      .recombWeights =     string choice of recombination weighting strategy
%*******************************************************************************
function guiCMAES = cmaes_gui(k,gui_var,gui)

% Set up "NW Options" figure window
guiCMAES.fig = figure(...
   'Name', ['CMAES Toolbox Options (Search #', int2str(k),')'], ...
   'Tag',  'CMAES_GUI', ...
   'DefaultUIControlUnits',     'normalized', ...
   'DefaultUIControlFontUnits', 'normalized', ...
   'DefaultUIControlFontName',  'Helvetica',  ...
   'DefaultUIControlFontSize',  0.5, ...
   'DefaultUIControlHorizontalAlignment', 'left', ...
   'DefaultUIControlStyle',     'text', ...
   'WindowStyle',               'modal', ...
   'Resize',                    'on', ...
   'Units',                     'normalized', ...
   'Position',                  [0.20 0.20 0.70 0.70], ...
   'MenuBar',                   'none', ...
   'NumberTitle',               'off');

% Figure window shading to make it appear recessed
uipanel;
uicontrol(guiCMAES.fig, 'Style','frame','Position', [0, 0, 1, .004],    ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(guiCMAES.fig, 'Style','frame','Position', [.997, 0, .003, 1], ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(guiCMAES.fig, 'Style','frame','Position', [0, .996, 1, .004], ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
uicontrol(guiCMAES.fig, 'Style','frame','Position', [0, 0, .003, 1],    ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);

% Labels for CMAES screen objects and data fields
uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              'Population Size: ', ...
    'Position',            [.07 .895 .15 .05]);
guiCMAES.popSize         = uicontrol(guiCMAES.fig, ...
    'Style',               'edit', ...
    'String',              gui_var.Options.cmaes(k).popSize, ...
    'BackgroundColor',     'white', ...
    'HorizontalAlignment', 'center', ...
    'Position',            [.22 .90 .20 .05]);
 
uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              'Population Size Increment: ', ...
    'Position',            [.52 .895 .20 .05]);
guiCMAES.incPopSize      = uicontrol(guiCMAES.fig, ...
    'Style',               'edit', ...
    'String',              gui_var.Options.cmaes(k).incPopSize, ...
    'BackgroundColor',     'white', ...
    'HorizontalAlignment', 'center', ...
    'Position',            [.72 .90 .20 .05]);

uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              '# Parents: ', ...
    'Position',            [.07 .805 .15 .05]);
guiCMAES.nParents        = uicontrol(guiCMAES.fig, ...
    'Style',               'edit', ...
    'String',              gui_var.Options.cmaes(k).nParents, ...
    'BackgroundColor',     'white', ...
    'HorizontalAlignment', 'center', ...
    'Position',            [.22 .81 .20 .05]);

uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              '# Restarts: ', ...
    'Position',            [.52 .805 .20 .05]);
guiCMAES.nRestarts       = uicontrol(guiCMAES.fig, ...
    'Style',               'edit', ...
    'String',              gui_var.Options.cmaes(k).nRestarts, ...
    'BackgroundColor',     'white', ...
    'HorizontalAlignment', 'center', ...
    'Position',            [.72 .81 .20 .05]);

uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              'Standard Deviations: ', ...
    'Position',            [.07 .715 .15 .05]);
guiCMAES.sigma           = uicontrol(guiCMAES.fig, ...
    'Style',               'edit', ...
    'String',              gui_var.Options.cmaes(k).sigma, ...
    'BackgroundColor',     'white', ...
    'HorizontalAlignment', 'center', ...
    'Position',            [.22 .72 .20 .05]);

uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              'Recombination Weighting: ', ...
    'Position',            [.52 .715 .20 .05]);
guiCMAES.recombWeights   = uicontrol(guiCMAES.fig, ...
    'Style',               'popupmenu', ...
    'String',              gui_var.Labels.cmaesRecombWeights, ...
    'BackgroundColor',     'white', ...
    'Value',               gui_var.Choice.cmaesRecombWeights(k), ...
    'Position',            [.72 .72 .20 .05]);

uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              'STOPPING CRITERIA', ...
    'HorizontalAlignment', 'center', ...
    'FontSize',            .5, ...
    'Position',            [.40 .60 .20 .05]);

uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              'Fitness function lower bound: ', ...
    'Position',            [.30 .525 .20 .05]);
guiCMAES.Stop.fitness    = uicontrol(guiCMAES.fig, ...
    'Style',               'edit', ...
    'String',              gui_var.Options.cmaes(k).Stop.fitness, ...
    'BackgroundColor',     'white', ...
    'HorizontalAlignment', 'center', ...
    'Position',            [.52 .53 .20 .05]);

uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              'Max # of Function Evaluations: ', ...
    'Position',            [.30 .455 .20 .05]);
guiCMAES.Stop.maxFunEval = uicontrol(guiCMAES.fig, ...
    'Style',               'edit', ...
    'String',              gui_var.Options.cmaes(k).Stop.maxFunEval, ...
    'BackgroundColor',     'white', ...
    'HorizontalAlignment', 'center', ...
    'Position',            [.52 .46 .20 .05]);

uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              'Max # of Iterations: ', ...
    'Position',            [.30 .385 .20 .05]);
guiCMAES.Stop.maxIter = uicontrol(guiCMAES.fig, ...
    'Style',               'edit', ...
    'String',              gui_var.Options.cmaes(k).Stop.maxIter, ...
    'BackgroundColor',     'white', ...
    'HorizontalAlignment', 'center', ...
    'Position',            [.52 .39 .20 .05]);

uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              'Minimal variable change: ', ...
    'Position',            [.30 .315 .20 .05]);
guiCMAES.Stop.tolX       = uicontrol(guiCMAES.fig, ...
    'Style',               'edit', ...
    'String',              gui_var.Options.cmaes(k).Stop.tolX, ...
    'BackgroundColor',     'white', ...
    'HorizontalAlignment', 'center', ...
    'Position',            [.52 .32 .20 .05]);

uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              'Maximal variable change: ', ...
    'Position',            [.30 .245 .20 .05]);
guiCMAES.Stop.tolUpX     = uicontrol(guiCMAES.fig, ...
    'Style',               'edit', ...
    'String',              gui_var.Options.cmaes(k).Stop.tolUpX, ...
    'BackgroundColor',     'white', ...
    'HorizontalAlignment', 'center', ...
    'Position',            [.52 .25 .20 .05]);

uicontrol(guiCMAES.fig, ...
    'Style',               'text', ...
    'String',              'Minimal objective change: ', ...
    'Position',            [.30 .175 .20 .05]);
guiCMAES.Stop.tolFun     = uicontrol(guiCMAES.fig, ...
    'Style',               'edit', ...
    'String',              gui_var.Options.cmaes(k).Stop.tolFun, ...
    'BackgroundColor',     'white', ...
    'HorizontalAlignment', 'center', ...
    'Position',            [.52 .18 .20 .05]);

% guiCMAES.merit             = uicontrol(guiCMAES.fig, ...
%       'Style',               'checkbox', ...
%       'String',              'Use Merit Function to Penalize Clustering', ...
%       'Value',               gui_var.Options.Search(k).merit, ...
%       'Position',            [.02 .21 .70 .09]);

% The Done and Cancel Buttons
%    'Position',   [.33, .04, .10, .08], ...
guiCMAES.done              = uicontrol(guiCMAES.fig, ...
      'Style',               'pushbutton', ...
      'String',              'Done', ...
      'FontWeight',          'bold', ...
      'Position',            [.30, .02, .12, .09], ...
      'Callback',            {gui.func.loadCMAESOptions,1,k});
guiCMAES.cancel            = uicontrol(guiCMAES.fig, ...
      'Style',               'pushbutton', ...
      'String',              'Cancel', ...
      'FontWeight',          'bold', ...
      'Position',            [.60, .02, .12, .09], ...
      'Callback',            {gui.func.loadCMAESOptions,0,k});
setappdata(guiCMAES.fig,'guiCMAES',guiCMAES);
end

%*******************************************************************************
% loadCMAESOptions:  Load CMAES Options.
% ------------------------------------------------------------------------------
% Called by: cmaes_gui
% VARIABLES:
%  loadCMAES           = flag indicating if CMAES options are accepted
%  k                   = Search number
%  figSearch           = handle for Search figure window
%  guiSearch           = structure of Search window object handles
%    .Screen           =   handles for Search figure window objects
%      .type           =     type of Search
%  figCMAES            = handle for CMAES figure window
%  guiCMAES            = structure of handles for CMAES window objects
%    .sigma            =   handle for CMAES sigma field
%    .Stop             =   handles for CMAES stopping criteria 
%      .fitness        =     handle for minimal fitness function value
%      .maxFunEvals    =     handle for maximum number of function evaluations
%      .maxIter        =     handle for maximum number of iterations
%      .tolX           =     handle for minimal variable change
%      .tolUpX         =     handle for maximal variable change
%      .tolFun         =     handle for minimal objective function value change
%    .evalInitialX     =   handle for checkbox to evaluate initial point
%    .incPopSize       =   handle for population size increment factor
%    .popSize          =   handle for initial population size
%    .nParents         =   handle for CMAES number of parents
%    .nRestarts        =   handle for number of restarts
%    .recombWeights    =   handle for recoombination weighting strategy
%  fig                 = handle for the NOMADm figure window
%  gui_var             = structure of all GUI variables
%    .Types            =   lists of possible types
%      .recombWeights  =     list of possible recombination weighting strategies
%    .Options.cmaes(k) =   structure of CMAES parameters for Search k
%  previousType        = previously selected Search type
%  cmaes               =   sturcture of user-chosen CMAES parameters
%    .sigma            =     standard deviations for each variable
%    .Stop             =     structure of stopping conditions
%      .Fitness        =       fitness value lower threshold
%      .MaxFunEvals    =       maximum number of function evaluations
%      .MaxIter        =       maximum number of iterations
%      .TolX           =       minimum distance between x-values
%      .TolUpX         =       maximum distance between x-values
%      .TolFun         =       minimum difference between function values
%    .evalInitialX     =     flag for evaluating the initial point
%    .nRestarts        =    number of restarts
%    .incPopSize       =     factor for multiplying population size
%    .popSize          =     population size 
%    .nParents         =     number of parents
%    .recombWeights    =     string choice of recombination weighting strategy
%*******************************************************************************
function loadCMAESOptions(~,~,loadCMAES,k) %#ok

figSearch = findobj('Tag','Search_GUI');
guiSearch = getappdata(figSearch,'guiSearch');
figCMAES  = findobj('Tag','CMAES_GUI');
guiCMAES  = getappdata(figCMAES,'guiCMAES');
fig       = findobj('Tag','NOMADm_GUI');
gui_var   = getappdata(fig,'gui_var');

if loadCMAES

   % Load default values first
   gui_var.Options.cmaes(k)   = gui_var.Defaults.Options.cmaes(k);

   % Load options from screen
   cmaes = gui_var.Options.cmaes(k);
   cmaes.sigma                = get(guiCMAES.sigma,           'String');
   cmaes.Stop.fitness         = get(guiCMAES.Stop.fitness,    'String');
   cmaes.Stop.maxFunEval      = get(guiCMAES.Stop.maxFunEval, 'String');
   cmaes.Stop.maxIter         = get(guiCMAES.Stop.maxIter,    'String');
   cmaes.Stop.tolX            = get(guiCMAES.Stop.tolX,       'String');
   cmaes.Stop.tolUpX          = get(guiCMAES.Stop.tolUpX,     'String');
   cmaes.Stop.tolFun          = get(guiCMAES.Stop.tolFun,     'String');
   cmaes.incPopSize           = get(guiCMAES.incPopSize,      'String');
   cmaes.popSize              = get(guiCMAES.popSize,         'String'); 
   cmaes.nParents             = get(guiCMAES.nParents,        'String');
   cmaes.nRestarts            = get(guiCMAES.nRestarts,       'String');
   cmaes.recombWeights        = gui_var.Types.cmaesRecombWeights{...
                                        get(guiCMAES.recombWeights,'Value')};
   gui_var.Options.cmaes(k)   = cmaes;
   setappdata(fig,'gui_var',gui_var);
else
   previousType = get(guiSearch.Screen(k).type,'UserData');
   set(guiSearch.Screen(k).type,'Value',previousType{2});
end
close(guiCMAES.fig);
end

%*******************************************************************************
% dace_gui:  Displays a user input screen to set DACE Toolbox parameters.
% ------------------------------------------------------------------------------
% Called by: loadUserSearchFile
% Calls:     loadSurrogateOptions
% VARIABLES:
%  k                    = Search number for this DACE screen
%  guiSurrogate         = structure of all DACE window object handles
%    .fig               =   DACE figure window
%    .daceRegr          =   regression function popup menu
%    .daceCorr          =   rcorrelation function popup menu
%    .daceIsotropic     =   checkbox for isotropic correlations
%    .maxCondR          =   field for max correlation matrix condition number
%    .upper             =   field for entering upper bound for theta
%    .local             =   checkbox for restricting model to local region
%    .merit             =   checkbox to use merit function
%    .done              =   the Done pushbutton
%    .cancel            =   the Cancel pushbutton
%  fig                  = handle for the NOMADm figure window
%  gui_var              = structure of all GUI variables
%    .Labels            =   labels for DACE function popup menus
%      .daceRegression  =     labels for DACE regression functions
%      .daceCorrelation =     labels for DACE correlation functions
%    .Choice            =   integer choices for DACE functions
%      .daceRegr        =     choice for DACE regression function
%      .daceCorr        =     choice for DACE correlation function
%  gui.func             = structure of all GUI function handles
%  Options              = structure of MADS options
%    .dace              =   user-chosen DACE Toolbox parameters
%      .isotropic       =   flag for isotropic theta
%*******************************************************************************
function guiSurrogate = dace_gui(k,gui_var,gui)

% Set up "DACE Options" figure window
guiSurrogate.fig = figure(...
   'Name', ['DACE Toolbox Options (Search #', int2str(k),')'], ...
   'Tag',  'Surrogate_GUI', ...
   'DefaultUIControlUnits',          'normalized', ...
   'DefaultUIControlFontUnits',      'normalized', ...
   'DefaultUIControlFontName',       'Helvetica',  ...
   'DefaultUIControlFontSize',        0.5, ...,
   'DefaultUIControlStyle',           'text', ...
   'WindowStyle',                     'modal', ...
   'Resize',                          'on', ...
   'Units',                           'normalized', ...
   'Position',                        [0.20 0.25 0.50 0.50], ...
   'MenuBar',                         'none', ...
   'NumberTitle',                     'off');

% Figure window shading to make it appear recessed
uipanel;
uicontrol(guiSurrogate.fig, 'Style','frame','Position', [0, 0, 1, .004],    ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(guiSurrogate.fig, 'Style','frame','Position', [.997, 0, .003, 1], ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(guiSurrogate.fig, 'Style','frame','Position', [0, .996, 1, .004], ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
uicontrol(guiSurrogate.fig, 'Style','frame','Position', [0, 0, .003, 1],    ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);

% Labels for DACE screen objects and data fields
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'Regression Model: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .885 .34 .06]);
guiSurrogate.daceRegr      = uicontrol(guiSurrogate.fig, ...
      'Style',               'popupmenu', ...
      'String',              gui_var.Labels.daceRegression, ...
      'BackgroundColor',     'white', ...
      'Value',               gui_var.Choice.daceRegr(k), ...
      'Position',            [.37 .89 .60 .06]);
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'Correlation Model: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .815 .34 .06]);
guiSurrogate.daceCorr      = uicontrol(guiSurrogate.fig, ...
      'Style',               'popupmenu', ...
      'String',              gui_var.Labels.daceCorrelation, ...
      'BackgroundColor',     'white', ...
      'Value',               gui_var.Choice.daceCorr(k), ...
      'Position',            [.37 .82 .60 .06]);
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'MLE Optimization Scheme: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .745 .34 .06]);
guiSurrogate.mleOptType    = uicontrol(guiSurrogate.fig, ...
      'Style',               'popupmenu', ...
      'String',              gui_var.Labels.mleOpt, ...
      'BackgroundColor',     'white', ...
      'Value',               gui_var.Choice.mleOpt(k), ...
      'Position',            [.37 .75 .60 .06]);
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'Maximum Correlation Matrix Condition Number: ',...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .675 .64 .06]);
guiSurrogate.maxCondR   = uicontrol(guiSurrogate.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.dace(k).maxCondR), ...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.67 .68 .30 .06]);
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'Minimum Log-Likelihood Derivative Norm: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .605 .64 .06]);
guiSurrogate.minNormDL     = uicontrol(guiSurrogate.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.dace(k).minNormDL), ...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.67 .61 .30 .06]);
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'MLE Bounds Tolerance: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .535 .64 .06]);
guiSurrogate.mleTolBounds  = uicontrol(guiSurrogate.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.dace(k).mleTolBounds), ...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.67 .54 .30 .06]);
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'MLE # Grid/MultiStart Points: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .465 .64 .06]);
guiSurrogate.mleOptPoints  = uicontrol(guiSurrogate.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.dace(k).mleOpt.nPoints),...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.67 .47 .30 .06]);
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'CBLGS Points: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .395 .64 .06]);
guiSurrogate.cblgsGoal     = uicontrol(guiSurrogate.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.Search(k).cblgs.nGoal),...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.67 .40 .30 .06]);
guiSurrogate.daceIsotropic = uicontrol(guiSurrogate.fig, ...
      'Style',               'checkbox',  ...
      'String',              'Isotropic Fitting Parameters', ...
      'Value',               gui_var.Options.dace(k).isotropic,   ...
      'Position',            [.02 .32 .70 .07]);
guiSurrogate.local         = uicontrol(guiSurrogate.fig, ...
      'Style',               'checkbox', ...
      'String',              'Restrict Model to Trust Region', ...
      'Value',               gui_var.Options.Search(k).local, ...
      'Position',            [.02 .25 .70 .07]);
guiSurrogate.merit              = uicontrol(guiSurrogate.fig, ...
      'Style',               'checkbox', ...
      'String',              'Use Merit Function to Penalize Clustering', ...
      'Value',               gui_var.Options.Search(k).merit, ...
      'Position',            [.02 .18 .70 .07]);

% The Done and Cancel Buttons
guiSurrogate.done          = uicontrol(guiSurrogate.fig, ...
      'Style',               'pushbutton', ...
      'String',              'Done', ...
      'FontWeight',          'bold', ...
      'Position',            [.27, .02, .18, .10], ...
      'Callback',            {gui.func.loadSurrogateOptions,'DACE',k});
guiSurrogate.cancel        = uicontrol(guiSurrogate.fig, ...
      'Style',               'pushbutton', ...
      'String',              'Cancel', ...
      'FontWeight',          'bold', ...
      'Position',            [.54, .02, .18, .10], ...
      'Callback',            {gui.func.loadSurrogateOptions,'',k});
setappdata(guiSurrogate.fig, 'guiSurrogate',guiSurrogate);
end

%*******************************************************************************
% nw_gui:  Displays a user input screen to set NW Toolbox parameters.
% ------------------------------------------------------------------------------
% Called by: loadUserSearchFile
% Calls:     loadSurrogateOptions
% VARIABLES:
%  k                    = Search number for this NW screen
%  guiSurrogate         = structure of all NW window object handles
%    .fig               =   NW figure window
%    .nwKernel          =   kernel function popup menu
%    .nwSigma           =   sigma parameter initial guess 
%    .lower             =   sigma parameter lower bound
%    .upper             =   sigma parameter upper bound
%    .local             =   checkbox for restricting model to local region
%    .merit             =   checkbox to use merit function
%    .done              =   the Done pushbutton
%    .cancel            =   the Cancel pushbutton
%  fig                  = handle for the NOMADm figure window
%  gui_var              = structure of all GUI variables
%    .Labels            =   labels for NW function popup menus
%      .nwKernel        =     labels for NW Kernel functions
%    .Choice            =   integer choices for NW functions
%      .nwKernel        =     choice for NW kernel function
%  gui.func             = structure of all GUI function handles
%  Options              = structure of MADS options
%    .nw                =   user-chosen NW Toolbox parameters
%      .kernel          =   kernel function
%      .sigma           =   sigma parameter
%      .lower           =   sigma parameter lower bound
%      .upper           =   sigma parameter upper bound
%*******************************************************************************
function guiSurrogate = nw_gui(k,gui_var,gui)

% Set up "NW Options" figure window
guiSurrogate.fig = figure(...
   'Name', ['NW Toolbox Options (Search #', int2str(k),')'], ...
   'Tag',  'Surrogate_GUI', ...
   'DefaultUIControlUnits',     'normalized', ...
   'DefaultUIControlFontUnits', 'normalized', ...
   'DefaultUIControlFontName',  'Helvetica',  ...
   'DefaultUIControlFontSize',   0.5, ...,
   'DefaultUIControlStyle',     'text', ...
   'WindowStyle',               'modal', ...
   'Resize',                    'on', ...
   'Units',                     'normalized', ...
   'Position',                   [0.25 0.30 0.40 0.40], ...
   'MenuBar',                   'none', ...
   'NumberTitle',               'off');

% Figure window shading to make it appear recessed
uipanel;
uicontrol(guiSurrogate.fig, 'Style','frame','Position', [0, 0, 1, .004],    ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(guiSurrogate.fig, 'Style','frame','Position', [.997, 0, .003, 1], ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(guiSurrogate.fig, 'Style','frame','Position', [0, .996, 1, .004], ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
uicontrol(guiSurrogate.fig, 'Style','frame','Position', [0, 0, .003, 1],    ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);

% Labels for NW screen objects and data fields
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'Kernel Function: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .855 .34 .08]);
guiSurrogate.nwKernel      = uicontrol(guiSurrogate.fig, ...
      'Style',               'popupmenu', ...
      'String',              gui_var.Labels.nwKernel, ...
      'BackgroundColor',     'white', ...
      'Value',               gui_var.Choice.nwKernel(k), ...
      'Position',            [.37 .86 .60 .09]);
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'Sigma Initial Guess: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .675 .34 .08]);
guiSurrogate.nwSigma       = uicontrol(guiSurrogate.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.nw(k).sigma), ...
      'FontSize',            .6, ...
      'BackgroundColor',     'white', ...
      'Position',            [.37 .69 .35 .08]);
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'Sigma Lower Bound: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .585 .34 .08]);
guiSurrogate.lower         = uicontrol(guiSurrogate.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.nw(k).lower), ...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.37 .60 .35 .08]);
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'Sigma Upper Bound: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .495 .34 .08]);
guiSurrogate.upper         = uicontrol(guiSurrogate.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.nw(k).upper), ...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.37 .51 .35 .08]);
guiSurrogate.local         = uicontrol(guiSurrogate.fig, ...
      'Style',               'checkbox', ...
      'String',              'Restrict Model to Trust Region', ...
      'Value',               gui_var.Options.Search(k).local, ...
      'Position',            [.02 .30 .70 .09]);
guiSurrogate.merit         = uicontrol(guiSurrogate.fig, ...
      'Style',               'checkbox', ...
      'String',              'Use Merit Function to Penalize Clustering', ...
      'Value',               gui_var.Options.Search(k).merit, ...
      'Position',            [.02 .21 .70 .09]);

% The Done and Cancel Buttons
guiSurrogate.done          = uicontrol(guiSurrogate.fig, ...
      'Style',               'pushbutton', ...
      'String',              'Done', ...
      'FontWeight',          'bold', ...
      'Position',            [.27, .02, .18, .12], ...
      'Callback',            {gui.func.loadSurrogateOptions,'NW',k});
guiSurrogate.cancel        = uicontrol(guiSurrogate.fig, ...
      'Style',               'pushbutton', ...
      'String',              'Cancel', ...
      'FontWeight',          'bold', ...
      'Position',            [.54, .02, .18, .12], ...
      'Callback',            {gui.func.loadSurrogateOptions,'',k});
setappdata(guiSurrogate.fig, 'guiSurrogate',guiSurrogate);
end

%*******************************************************************************
% rbf_gui:  Displays a user input screen to set NW Toolbox parameters.
% ------------------------------------------------------------------------------
% Called by: loadUserSearchFile
% Calls:     loadSurrogateOptions
% VARIABLES:
%  k                    = Search number for this RBF screen
%  guiSurrogate         = structure of all RBF window object handles
%    .fig               =   RBF figure window
%    .rbfKernel         =   kernel function popup menu
%    .local             =   checkbox for restricting model to local region
%    .merit             =   checkbox to use merit function
%    .done              =   the Done pushbutton
%    .cancel            =   the Cancel pushbutton
%  fig                  = handle for the NOMADm figure window
%  gui_var              = structure of all GUI variables
%    .Labels            =   labels for RBF function popup menus
%      .rbfKernel       =     labels for RBF Kernel functions
%    .Choice            =   integer choices for RBF functions
%      .rbfKernel       =     choice for RBF kernel function
%  gui.func             = structure of all GUI function handles
%  Options              = structure of MADS options
%    .nw                =   user-chosen RBF Toolbox parameters
%      .kernel          =   kernel function
%*******************************************************************************
function guiSurrogate = rbf_gui(k,gui_var,gui)

% Set up "NW Options" figure window
guiSurrogate.fig = figure(...
   'Name', ['RBF Toolbox Options (Search #', int2str(k),')'], ...
   'Tag',  'Surrogate_GUI', ...
   'DefaultUIControlUnits',          'normalized', ...
   'DefaultUIControlFontUnits',      'normalized', ...
   'DefaultUIControlFontName',       'Helvetica',  ...
   'DefaultUIControlFontSize',        0.5, ...,
   'DefaultUIControlStyle',           'text', ...
   'WindowStyle',                     'modal', ...
   'Resize',                          'on', ...
   'Units',                           'normalized', ...
   'Position',                        [0.25 0.30 0.40 0.40], ...
   'MenuBar',                         'none', ...
   'NumberTitle',                     'off');

% Figure window shading to make it appear recessed
uipanel;
uicontrol(guiSurrogate.fig, 'Style','frame','Position', [0, 0, 1, .004],    ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(guiSurrogate.fig, 'Style','frame','Position', [.997, 0, .003, 1], ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(guiSurrogate.fig, 'Style','frame','Position', [0, .996, 1, .004], ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
uicontrol(guiSurrogate.fig, 'Style','frame','Position', [0, 0, .003, 1],    ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);

% Labels for RBF screen objects and data fields
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'Kernel Function: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .855 .34 .08]);
guiSurrogate.rbfKernel     = uicontrol(guiSurrogate.fig, ...
      'Style',               'popupmenu', ...
      'String',              gui_var.Labels.rbfKernel, ...
      'BackgroundColor',     'white', ...
      'Value',               gui_var.Choice.rbfKernel(k), ...
      'Position',            [.37 .86 .60 .09]);
uicontrol(guiSurrogate.fig, ...
      'Style',               'text', ...
      'String',              'Polynomial Function: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .735 .34 .08]);
guiSurrogate.rbfPoly       = uicontrol(guiSurrogate.fig, ...
      'Style',               'popupmenu', ...
      'String',              gui_var.Labels.rbfPoly, ...
      'BackgroundColor',     'white', ...
      'Value',               gui_var.Choice.rbfPoly(k), ...
      'Position',            [.37 .74 .60 .09]);
guiSurrogate.local         = uicontrol(guiSurrogate.fig, ...
      'Style',               'checkbox', ...
      'String',              'Restrict Model to Trust Region', ...
      'Value',               gui_var.Options.Search(k).local, ...
      'Position',            [.02 .30 .70 .09]);
guiSurrogate.merit         = uicontrol(guiSurrogate.fig, ...
      'Style',               'checkbox', ...
      'String',              'Use Merit Function to Penalize Clustering', ...
      'Value',               gui_var.Options.Search(k).merit, ...
      'Position',            [.02 .21 .70 .09]);

% The Done and Cancel Buttons
guiSurrogate.done          = uicontrol(guiSurrogate.fig, ...
      'Style',               'pushbutton', ...
      'String',              'Done', ...
      'FontWeight',          'bold', ...
      'Position',            [.27, .02, .18, .12], ...
      'Callback',            {gui.func.loadSurrogateOptions,'RBF',k});
guiSurrogate.cancel        = uicontrol(guiSurrogate.fig, ...
      'Style',               'pushbutton', ...
      'String',              'Cancel', ...
      'FontWeight',          'bold', ...
      'Position',            [.54, .02, .18, .12], ...
      'Callback',            {gui.func.loadSurrogateOptions,'',k});
setappdata(guiSurrogate.fig, 'guiSurrogate',guiSurrogate);
end

%*******************************************************************************
% loadSurrogateOptions: Select and load surrogate options.
% ------------------------------------------------------------------------------
% Called by: dace_gui, nw_gui, rbf_gui (callback: Done and Cancel pushbuttons)
% Calls:     modifySearchScreen
% VARIABLES:
%  loadSurrogate      = flag indicating if Surrogate options are accepted
%  k                  = Search number
%  figSearch          = handle for Search figure window
%  guiSearch          = structure of Search window object handles
%    .Screen          =   handles for Search figure window objects
%      .type          =     type of Search
%      .sfile         =     names of optional surrogate files
%  figSurrogate       = handle for Surrogate figure window
%  guiSurrogate       = structure of handles for Surrogate window objects
%    .daceRegr        =   handle for DACE regression function popup menu
%    .daceCorr        =   handle for DACE correlation function popup menu
%    .daceTheta       =   handle for DACE correlation parameter initial guess
%    .daceIsotropic   =   handle for DACE isotropic parameters check box
%    .daceSetBounds   =   handle for DACE theta bound manual setting check box
%    .nwKernel        =   handle for NW kernel function popup menu
%    .nwSigma         =   handle for NW sigma parameter initial guess 
%    .lower           =   handle for surrogate parameter lower bound
%    .upper           =   handle for surrogate parameter upper bound
%    .rbfKernel       =   handle for RBF kernel function popup menu
%    .local           =   checkbox handle for using a trust region
%    .merit           =   checkbox handle for using a merit function
%  fig                = handle for the NOMADm figure window
%  gui_var            = structure of all GUI variables
%    .Types           =   lists of possible types
%      .daceRegr      =     list of possible DACE regression functions
%      .daceCorr      =     list of possible DACE correlation functions
%      .nwKernel      =     list of possible NW regression functions
%      .rbfKernel     =     list of possible RBF regression functions
%    .Options.dace(k) =   structure of DACE parameters for Search k
%    .Options.nw(k)   =   structure of NW parameters for Search k
%  previousType       = previously selected Search type
%  dace               = structure of DACE parameters
%    .reg             =   regression function handle
%    .corr            =   correlation function handle
%    .theta           =   initial guess for correlation parameters
%    .lower           =   lower bounds for correlation parameters
%    .upper           =   upper bounds for correlation parameters
%    .isotropic       =   flag for isotropic correlation parameters
%    .setBounds       =   flag for manually setting bounds on theta
%  nw                 = structure of NW parameters
%    .kernel          =   kernel function handle
%    .sigma           =   initial guess for sigma parameter
%    .lower           =   lower bounds for sigma parameter
%    .upper           =   upper bounds for sigma parameter
%  rbf                = structure of RBF parameters
%    .kernel          =   kernel function handle
%*******************************************************************************
function loadSurrogateOptions(~,~,loadSurrogate,k) %#ok

figSearch    = findobj('Tag','Search_GUI');
guiSearch    = getappdata(figSearch,'guiSearch');
figSurrogate = findobj('Tag','Surrogate_GUI');
guiSurrogate = getappdata(figSurrogate,'guiSurrogate');
fig          = findobj('Tag','NOMADm_GUI');
gui_var      = getappdata(fig,'gui_var');

if isempty(loadSurrogate)
   previousType = get(guiSearch.Screen(k).type,'UserData');
   set(guiSearch.Screen(k).type,'Value',previousType{2});
else

   % Load default values first
   gui_var.Options.dace(k) = gui_var.Defaults.Options.dace(k);
   gui_var.Options.nw(k)   = gui_var.Defaults.Options.nw(k);
   gui_var.Options.rbf(k)  = gui_var.Defaults.Options.rbf(k);

   % Load surrogate-specific options
   switch upper(loadSurrogate)
   case 'DACE'
      dace        = gui_var.Defaults.Options.dace(k);
      dace.reg    = gui_var.Types.daceRegr{get(guiSurrogate.daceRegr,'Value')};
      dace.corr   = gui_var.Types.daceCorr{get(guiSurrogate.daceCorr,'Value')};
      dace.mleOpt.type = gui_var.Types.mleOpt{get(guiSurrogate.mleOptType,'Value')};
      dace.mleOpt.nPoints = str2double(get(guiSurrogate.mleOptPoints,'String'));
      dace.maxCondR       = str2double(get(guiSurrogate.maxCondR,    'String'));
      dace.minNormDL      = str2double(get(guiSurrogate.minNormDL,   'String'));
      dace.mleTolBounds   = str2double(get(guiSurrogate.mleTolBounds,'String'));
%       dace.cblgs.nGoal    = str2double(get(guiSurrogate.cblgsGoal,   'String'));
      dace.isotropic      = get(guiSurrogate.daceIsotropic,'Value');
      gui_var.Options.dace(k) = dace;
      gui_var.Options.Search(k).cblgs.nGoal = ...
              str2double(get(guiSurrogate.cblgsGoal,'String'));
      set(guiSearch.Screen(k).sfile,'UserData','', ...
                                    'String',[dace.reg,'/',dace.corr]);
   case 'NW'
      nw.kernel = gui_var.Types.nwKernel{get(guiSurrogate.nwKernel,'Value')};
      nw.sigma  = str2double(get(guiSurrogate.nwSigma,'String'));
      nw.lower  = str2double(get(guiSurrogate.lower,'String'));
      nw.upper  = str2double(get(guiSurrogate.upper,'String'));
      gui_var.Options.nw(k) = nw;
      set(guiSearch.Screen(k).sfile,'UserData','','String',nw.kernel);
   case 'RBF'
      rbf.kernel = gui_var.Types.rbfKernel{get(guiSurrogate.rbfKernel,'Value')};
      rbf.poly   = gui_var.Types.rbfPoly{get(guiSurrogate.rbfPoly,'Value')};
      gui_var.Options.rbf(k)  = rbf;
      set(guiSearch.Screen(k).sfile,'UserData','', ...
                                    'String',[rbf.kernel,'/',rbf.poly]);
   end
   setappdata(fig,'gui_var',gui_var);
   
   % Update certain Search screen fields and user data
   set(guiSearch.Screen(k).type, 'UserData', ...
                  {k,get(guiSearch.Screen(k).type,'Value'), ...
                     get(guiSurrogate.local,'Value'), ...
                     get(guiSurrogate.merit,'Value')});
   UserData = get(guiSurrogate.fig,'UserData');
   if isempty(UserData)
      set(guiSearch.Screen(k).file,'UserData','');
      set(guiSearch.Screen(k).file,'String',  '');
   else
      set(guiSearch.Screen(k).file,'UserData',UserData{1});
      set(guiSearch.Screen(k).file,'String',  UserData{2});
   end
end

% Close the surrogate window and update the Search screen
modifySearchScreen(2,0);
close(guiSurrogate.fig);
end

%*******************************************************************************
% updateScreen:  Displays a user input screen to set Search parameters.
% ------------------------------------------------------------------------------
% Called by: loadSession, resetSession
% Calls:     updateSearchLabels
% VARIABLES:
%  Choice               = structure of menu choices
%    .pollStrategy      =   integer choice of Poll strategy
%    .pollOrder         =   integer choice of Poll order strategy
%    .pollCenter        =   integer choice of Poll center
%  Options              = structure of MADS parameters
%    .nSearches         =   number of Search types used
%    .Search(n)         =   structure of Search parameters
%    .delta0            =   initial poll size
%    .deltaMax          =   maximum poll size
%    .meshRefine        =   mesh refinement factor
%    .meshCoarsen       =   mesh coarsening factor
%    .tolCache          =   tolerance for flagging point as being in Cache
%    .hmin              =   minimum h-value of an infeasible point
%    .hmax              =   maximum h-value of a filter point
%    .EPoll.fTrigger    =   f-value Extended Poll trigger
%    .EPoll.hTrigger    =   h-value Extended Poll trigger
%    .Term              =   substructure containing MADS termination criteria
%      .delta           =     poll size parameter
%      .iter            =     maximum number of MADS iterations
%      .func            =     maximum number of function evaluations
%      .time            =     maximum CPU time
%      .fails           =     maximum number of consecutive Poll failures
%    .TermFlag          =   substructure of termination criteria on/off switches
%      .iter            =     turns on/off number of iterations
%      .nFunc           =     turns on/off number of function evaluations
%      .time            =     turns on/off CPU time
%      .nFails          =     turns on/off number of consecutive Poll failures
%    .scale             =   flag for scaling mesh directions
%    .pollComplete      =   turns on/off complete Polling
%    .NPollComplete     =   turns on/off complete discrete neighbor Polling
%    .EPollComplete     =   turns on/off complete Extended Polling
%    .TermFlag.relative =   computes termination delta relative to .delta0
%    .useFilter         =   use filter for nonlinear constraints
%    .removeRedundancy  =   remove redundant linear constraints
%    .runStochastic     =   run as stochastic optimization problem
%    .fixCategorical    =   fix categorical variables as constant
%    .accelerate        =   flag for accelerating mesh refinement
%    .saveHistory       =   flag for saving history to a text file
%    .plotHistory       =   turns on/off a history plot
%    .plotFilter        =   turns on/off a real-time filter plot
%    .loadCache         =   flag for loading a pre-existing Cache file
%    .countCache        =   flag for counting Cache points as function calls
%    .runOneIteration   =   flag for running one MADS iteration at a time
%    .runUntilFeasible  =   flag for running MADS only until feasible
%  gui_var.Labels       = sub-structure of GUI labels
%    .pollStrategy      =   Poll strategy label
%    .pollCenter        =   Poll center label
%    .pollOrder         =   Poll order label
%  gui                  = structure of GUI object handles
%    .SearchType(k)     =   current Search types
%    .pollStrategy      =   current Poll strategy
%    .pollCenter        =   current Poll center
%    .pollOrder         =   current Poll order type
%    .delta0            =   current initial poll size
%    .deltaMax          =   current maximum poll size
%    .meshRefine        =   current mesh refinement factor
%    .meshCoarsen       =   current mesh coarsening factor
%    .tolCache          =   current Cache tolerance
%    .hmin              =   current minimum infeasible h-value
%    .hmax              =   current maximum filter h-value
%    .ePollXiF          =   current f-value Extended Poll trigger
%    .ePollXiH          =   current h-value Extended Poll trigger
%    .TermDelta         =   current mesh size termination criteria
%    .TermIter          =   current maximum number of MADS iterations
%    .TermFunc          =   current maximum number of function evaluations
%    .TermTime          =   current maximum CPU time
%    .TermFails         =   current max number of consec Poll failures
%    .TermFlag          =   structure of handles for termination checkboxes
%      .nIter           =     handle for number of iterations
%      .nFunc           =     handle for number of function evaluations
%      .time            =     handle for CPU time
%      .nFails          =     handle for number of consecutive Poll failures
%    .CacheMenu         =   Cache menu items
%    .MADSMenu          =   MADS menu items
%    .OptionsMenu       =   Options menu items
%    .RunMenu           =   Run menu items
%  onoff                = cell array of two strings "on" and "off"
%  maxDisplay           = number of Search Types displayed on the main GUI
%  searchLabel          = Search label that appears on the main GUI
%*******************************************************************************
function updateScreen(Choice,Options)

fig     = findobj('Tag','NOMADm_GUI');
gui_var = getappdata(fig,'gui_var');
gui     = getappdata(fig,'gui');
onoff   = {'on','off'};

% Update GUI Search fields
maxDisplay  = length(gui.searchLabel);
searchLabel = updateSearchLabels(maxDisplay,Options.nSearches,Options.Search);
for k = 1:maxDisplay
   set(gui.searchLabel(k), 'String', searchLabel{k});
end

% Update GUI text objects and checkboxes
set(gui.pollStrategy,'String',gui_var.Labels.poll(Choice.Poll.type));
set(gui.pollCenter,  'String',gui_var.Labels.pollCenter(Choice.Poll.center));
set(gui.pollOrder,   'String',gui_var.Labels.pollOrder(Choice.Poll.order));
set(gui.delta0,      'String', num2str(Options.delta0,        '%5.5g'));
set(gui.deltaMax,    'String', num2str(Options.deltaMax,      '%5.3g'));
set(gui.meshRefine,  'String', num2str(Options.meshRefine,    '%1.3g'));
set(gui.meshCoarsen, 'String', num2str(Options.meshCoarsen,   '%3.1f'));
set(gui.tolCache,    'String', num2str(Options.tolCache,      '%1.6g'));
set(gui.hmin,        'String', num2str(Options.hmin,          '%1.5g'));
set(gui.hmax,        'String', num2str(Options.hmax,          '%1.2f'));
set(gui.hRho,        'String', num2str(Options.hRho,          '%3.5g'));
set(gui.nRuns,       'String', num2str(Options.nRuns,         '%3d'));
set(gui.ePollXiF,    'String', num2str(Options.EPoll.fTrigger,'%3.5g'));
set(gui.ePollXiH,    'String', num2str(Options.EPoll.hTrigger,'%3.5g'));
set(gui.Term.delta,  'String', num2str(Options.Term.delta,    '%3.0g'));
field = fieldnames(gui.Term);
for k = 1:length(field)
   set(gui.Term.(field{k}), 'String', Options.Term.(field{k}));
end
field = fieldnames(gui.TermFlag);
for k = 1:length(field)
   set(gui.TermFlag.(field{k}), 'Value', Options.TermFlag.(field{k}));
end

% Update menu checkmarks
for k = 1:length(gui.degeneracyMenu)
   set(gui.degeneracyMenu(k), 'Checked', onoff{1+~strcmp(...
      Options.degeneracyScheme, gui_var.Types.degeneracyScheme{k})});
end

set(gui.scaleMenu(1),      'Checked',onoff{1+~(Options.scale ==  0)});
set(gui.scaleMenu(2),      'Checked',onoff{1+~(Options.scale ==  2)});
set(gui.scaleMenu(3),      'Checked',onoff{1+~(Options.scale == 10)});
set(gui.filterMenu(1),     'Checked',onoff{1+~(Options.useFilter == 0)});
set(gui.filterMenu(2),     'Checked',onoff{1+~(Options.useFilter == 1)});
set(gui.filterMenu(3),     'Checked',onoff{1+~(Options.useFilter == 2)});
set(gui.multiFiMenu(1),    'Checked',onoff{1+~(Options.useMultiFi == 0)});
set(gui.multiFiMenu(2),    'Checked',onoff{1+~(Options.useMultiFi == 1)});
set(gui.multiFiMenu(3),    'Checked',onoff{1+~(Options.useMultiFi == 2)});
set(gui.sensorMenu(1),     'Checked',onoff{1+~(Options.dimSensor == 0)});
set(gui.sensorMenu(2),     'Checked',onoff{1+~(Options.dimSensor == 1)});
set(gui.sensorMenu(3),     'Checked',onoff{1+~(Options.dimSensor == 2)});
set(gui.sensorMenu(4),     'Checked',onoff{1+~(Options.dimSensor == 3)});
set(gui.debugMenu(1),      'Checked',onoff{1+~(Options.debug == 0)});
set(gui.debugMenu(2),      'Checked',onoff{1+~(Options.debug == 1)});
set(gui.debugMenu(3),      'Checked',onoff{1+~(Options.debug == 2)});
set(gui.debugMenu(4),      'Checked',onoff{1+~(Options.debug == 3)});
set(gui.plotHistoryMenu(1),'Checked',onoff{1+~(Options.plotHistory == 0)});
set(gui.plotHistoryMenu(2),'Checked',onoff{1+~(Options.plotHistory == 1)});
set(gui.plotHistoryMenu(3),'Checked',onoff{1+~(Options.plotHistory == 2)});
set(gui.CacheMenu.load,            'Checked',onoff{2-Options.loadCache});
set(gui.CacheMenu.count,           'Checked',onoff{2-Options.countCache});
set(gui.CacheMenu.countNotInX,     'Checked',onoff{2-Options.countNotInX});
set(gui.OptionsMenu.runStochastic, 'Checked',onoff{2-Options.runStochastic});
set(gui.OptionsMenu.fixCategorical,'Checked',onoff{2-Options.fixCategorical});
set(gui.OptionsMenu.accelerate,    'Checked',onoff{2-Options.accelerate});
set(gui.OptionsMenu.saveHistory,   'Checked',onoff{2-Options.saveHistory});
set(gui.OptionsMenu.plotFilter,    'Checked',onoff{2-Options.plotFilter});
set(gui.RunMenu.oneIteration,      'Checked',onoff{2-Options.runOneIteration});
set(gui.RunMenu.execFeasible,      'Checked',onoff{2-Options.runUntilFeasible});
set(gui.OptionsMenu.TermFlag.relative, 'Checked', ...
    onoff{2-Options.TermFlag.relative});
set(gui.OptionsMenu.removeRedundancy,  'Checked', ...
    onoff{2-Options.removeRedundancy});
end

%*******************************************************************************
% updateSearchLabels:  Displays the Search Labels on the main GUI.
% ------------------------------------------------------------------------------
% Called by: loadSearchOptions, updateScreen
% VARIABLES:
%  label          = cell array containing labels for each Search type
%  maxDisplay     = number of Search types shown on the NOMADm GUI
%  nSearches      = number of Search types
%  Search(k)      = structure that describes the k-th Search
%    .type        =   string identifying Search type
%    .label       =   long text label for Search type
%    .nIter       =   number of iterations to perform Search
%    .nPoints     =   number of Search points
%    .sfile       =   names of optional user surrogate files
%    .file        =   name of  optional user file defining Search
%  k              = Search type counter
%  searchFilePath = path of optional Search file
%  searchFile     = name of optional Search file
%*******************************************************************************
function label = updateSearchLabels(maxDisplay,nSearches,Search)

% Initialize display parameters
[label{1:maxDisplay}] = deal('None');
if (nSearches > maxDisplay)
   label{maxDisplay} = 'Multiple Search Types';
   nSearches = maxDisplay - 1;
end

% Construct Search labels
for k = 1:nSearches
   if isfield(Search(k),'sfile')
      [~,sFile1] = fileparts(Search(k).sfile{1});
      [~,sFile2] = fileparts(Search(k).sfile{2});
   end
%   [~,uFile] = fileparts(Search(k).file);
   if ischar(Search(k).file)
      uFile = '';
   else
      uFile = func2str(Search(k).file);
   end
   switch Search(k).type
   case {'None','CCD'}
      label{k} = Search(k).label;
   case {'LHS','Mesh','CMAES','PS'}
      label{k} = [int2str(Search(k).nPoints),'-pt. ',Search(k).label];
   case {'SPollI','GPollI'}
      if (Search(k).nPoints == 1)
         label{k} = strrep(Search(k).label(1:end-1), 'N ', '');
      else
         label{k} = strrep(Search(k).label,'N',int2str(Search(k).nPoints));
      end
   case {'DACE'}
      label{k} = ['DACE Surrogate: ',   sFile1 ' / ', sFile2];
   case {'NW'}
      label{k} = ['N-W Surrogate: ',    sFile1 ' / ', sFile2];
   case {'RBF'}
      label{k} = ['RBF Surrogate: ',    sFile1 ' / ', sFile2];
   case {'SPS-DACE'}
      label{k} = ['SPS / DACE: ',       uFile, ' / ', sFile1 ' / ', sFile2];
   case {'SPS-NW'}
      label{k} = ['SPS / N-W: ',        uFile, ' / ', sFile1 ' / ', sFile2];
   case {'SPS-RBF'}
      label{k} = ['SPS / RBF: ',        uFile, ' / ', sFile1 ' / ', sFile2];
   case {'Custom'}
      label{k} = ['Custom Search: ',    uFile];
   case {'CustomS'}
      label{k} = ['Custom Surrogate: ', uFile];
   end
   label{k} = [label{k}, ' (', int2str(Search(k).nIter), ')'];
end
end

%*******************************************************************************
% results_gui:  Displays Results from MADS run.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Results-->View Run #k)
% Calls:     solutions_gui
% VARIABLES:
%  k               = run number
%  gui_var         = structure of GUI variables
%    .digits       =   number of digits accuracy to present on screen
%  gui.func        = structure of NOMADm figure window function handles
%  gui             = structure of NOMADm window object handles 
%    .problem      =   display of the optimization problem name
%    .Term         =   handles of termination criteria display
%      .delta      =     display of poll size termination criteria
%      .nIter      =     display of maximum number of MADS iterations
%      .nFunc      =     display of maximum number of function evaluations
%      .time       =     display of maximum CPU time
%      .nFails     =     display of maximum number of consecutive Poll failures
%  RunSet          = structure containing all information about a MADS run
%    .BestF        =   best feasible iterate for Run k
%    .BestI        =   least infeasible iterate for Run k
%    .RunStats     =   MADS Run k statistics (iterations, CPU time, etc.)
%      .delta      =     final poll size
%      .nIter      =     total number of iterations
%      .nFunc      =     total number of function evaluations
%      .time       =     total CPU time expended
%      .nFails     =     final number of consecutive Poll failures
%      .stopRun    =     final value of Stop Run button
%      .nCacheHits =     total number of Cache hits
%  k               = run number
%  noYes           = cell array containing the strings "No" and "Yes"
%  BestF, BestI    = best feasible and least infeasible solutions found
%  fig             = handle for the Results figure window
%  Run             = structure of handles for run statistics display
%    .delta        =   final poll size
%    .nIter        =   number of iterations
%    .nFunc        =   number of function evaluations
%    .nGrad        =   number of gradient evaluations
%    .time         =   CPU time expended
%    .nFails       =   number of consecutive Poll failures
%    .stop         =   interrupted by user (Yes/No)
%    .nCacheHits   =   number of Cache hits
%*******************************************************************************
function results_gui(~,~) %#ok

% Retrieve application data and initialize
gui_var = getappdata(gcbf,'gui_var');
gui     = getappdata(gcbf,'gui');
RunSet  = getappdata(gcbf,'RunSet');
k       = get(gcbo,'Position');
noYes   = {'No','Yes'};
BestF   = RunSet(k).BestF;
BestI   = RunSet(k).BestI;

% Set up "Display Results" figure window
fig = figure( ...
   'Name',                      ['MADS Results for Problem: ', ...
                                get(gui.problem, 'String'), ...
                                ', Run #' int2str(k)], ...
   'DefaultUIControlUnits',     'normalized', ...
   'DefaultUIControlFontUnits', 'normalized', ...
   'DefaultUIControlFontName',  'Helvetica',  ...
   'DefaultUIControlFontSize',   0.375, ...
   'DefaultUIControlStyle',      'text', ...
   'Units',                      'normalized', ...
   'Position',                   [0.1 0.15 0.8 0.7], ...
   'MenuBar',                    'none', ...
   'NumberTitle',                'off');

% Set up panels on the figure window
uipanel('Parent',fig,'Position',[0,  0,  1,.68],'BorderType','beveledin');
uipanel('Parent',fig,'Position',[0, .68,.5,.37],'BorderType','beveledin');
uipanel('Parent',fig,'Position',[.5,.68,.5,.37],'BorderType','beveledin');

% Display Best Feasible Solution with pushbutton to view the solution
uicontrol(fig, 'String','BEST FEASIBLE SOLUTION', ...
         'FontWeight','bold', ...
         'HorizontalAlignment','center','Position',[.02 .88 .46 .10]);

if (isempty(RunSet(k).BestF))
   uicontrol(fig, 'String','NO FEASIBLE ITERATES FOUND', ...
         'HorizontalAlignment','center','Position',[.02 .78 .46 .08]);
else
   uicontrol(fig, 'String','Objective Function Value: ', ...
         'HorizontalAlignment','left',  'Position',[.02 .82 .28 .08]);
   uicontrol(fig, 'String','Constraint Violation Measure: ', ...
         'HorizontalAlignment','left',  'Position',[.02 .77 .28 .08]);
   uicontrol(fig, 'String',num2str(RunSet(k).BestF.f,gui_var.Options.digits),...
         'HorizontalAlignment','right', 'Position',[.30 .82 .18 .08]);
   uicontrol(fig, 'String',num2str(RunSet(k).BestF.h,gui_var.Options.digits),...
         'HorizontalAlignment','right', 'Position',[.30 .77 .18 .08]);
   uicontrol(fig, 'String',             'View Solution', ...
                  'Style',              'pushbutton', ...
                  'FontWeight',         'bold', ...
                  'HorizontalAlignment','center', ...
                  'Position',           [.02 .71 .46 .08], ...
                  'UserData',           'Best Feasible', ...
                  'Callback',           {gui.func.solution_gui,k,BestF});
end

% Display Least Infeasible Solution with pushbutton to view the solution
if ~isempty(RunSet(k).BestI)

   uicontrol(fig, 'String','LEAST INFEASIBLE SOLUTION', ...
            'FontWeight','bold', ...
            'HorizontalAlignment','center','Position',[.52 .88 .46 .10]);
   uicontrol(fig,  'String','Objective Function Value: ', ...
         'HorizontalAlignment','left',  'Position',[.52 .82 .28 .08]);
   uicontrol(fig,'String','Constraint Violation Measure: ', ...
         'HorizontalAlignment','left',  'Position',[.52 .77 .28 .08]);
   uicontrol(fig, 'String',num2str(RunSet(k).BestI.f,gui_var.Options.digits),...
         'HorizontalAlignment','right', 'Position',[.80 .82 .18 .08]);
   uicontrol(fig, 'String',num2str(RunSet(k).BestI.h,gui_var.Options.digits),...
         'HorizontalAlignment','right', 'Position',[.80 .77 .18 .08]);
   uicontrol(fig, 'String',             'View Solution', ...
                  'Style',              'pushbutton', ...
                  'FontWeight',         'bold', ...
                  'HorizontalAlignment','center', ...
                  'Position',           [.52 .71 .46 .08], ...
                  'UserData',           'Least Infeasible', ...
                  'Callback',           {gui.func.solution_gui,k,BestI});
end

% Display Labels for MADS Run Statistics
uicontrol(fig, 'String','RUN STATISTICS', ...
   'FontWeight','bold', ...
   'HorizontalAlignment','center','Position',[.30 .56 .40 .10]);
uicontrol(fig, 'String','Final Poll Size:', ...
   'HorizontalAlignment','left',  'Position',[.20 .50 .35 .08]);
uicontrol(fig, 'String','MADS Iterations:', ...
   'HorizontalAlignment','left',  'Position',[.20 .45 .35 .08]);
uicontrol(fig, 'String','Poll Steps Executed:', ...
   'HorizontalAlignment','left',  'Position',[.20 .40 .35 .08]);
uicontrol(fig, 'String','Consecutive Poll Failures:', ...
   'HorizontalAlignment','left',  'Position',[.20 .35 .35 .08]);
uicontrol(fig, 'String','Function Evaluations:', ...
   'HorizontalAlignment','left',  'Position',[.20 .30 .35 .08]);
uicontrol(fig, 'String','Gradient Evaluations:', ...
   'HorizontalAlignment','left',  'Position',[.20 .25 .35 .08]);
uicontrol(fig, 'String','CPU Time:', ...
   'HorizontalAlignment','left',  'Position',[.20 .20 .35 .08]);
uicontrol(fig, 'String','Cache Hits:', ...
   'HorizontalAlignment','left',  'Position',[.20 .15 .35 .08]);
uicontrol(fig, 'String','Interrupted by User:', ...
   'HorizontalAlignment','left',  'Position',[.20 .10 .35 .08]);

% Display MADS Run Statistics (Compare with termination criteria)
Run.delta       = uicontrol(fig, ...
   'String',num2str(RunSet(k).RunStats.delta,gui_var.Options.digits), ...
   'HorizontalAlignment','right',  'Position',[.55 .50 .20 .08]);
Run.nIter       = uicontrol(fig, ...
   'String',int2str(RunSet(k).RunStats.nIter), ...
   'HorizontalAlignment','right',  'Position',[.55 .45 .20 .08]);
Run.nPoll       = uicontrol(fig, ...
   'String',int2str(RunSet(k).RunStats.nPoll), ...
   'HorizontalAlignment','right',  'Position',[.55 .40 .20 .08]);
Run.nFails      = uicontrol(fig, ...
   'String',num2str(RunSet(k).RunStats.nFails), ...
   'HorizontalAlignment','right',  'Position',[.55 .35 .20 .08]);
Run.nFunc       = uicontrol(fig, ...
   'String',int2str(RunSet(k).RunStats.nFunc), ...
   'HorizontalAlignment','right',  'Position',[.55 .30 .20 .08]);
Run.nGrad       = uicontrol(fig, ...
   'String',int2str(RunSet(k).RunStats.grad), ...
   'HorizontalAlignment','right',  'Position',[.55 .25 .20 .08]);
Run.time        = uicontrol(fig, ...
   'String',num2str(RunSet(k).RunStats.time), ...
   'HorizontalAlignment','right',  'Position',[.55 .20 .20 .08]);
Run.nCacheHits  = uicontrol(fig, ...
   'String',num2str(RunSet(k).RunStats.nCacheHits), ...
   'HorizontalAlignment','right',  'Position',[.55 .15 .20 .08]);
Run.stop        = uicontrol(fig, ...
   'String',noYes{1+RunSet(k).RunStats.stopRun}, ...
   'HorizontalAlignment','right',  'Position',[.55 .10 .20 .08]);

% Change colors for violated Termination Criteria
if (RunSet(k).RunStats.delta < str2double(get(gui.Term.delta,'String')))
   set(Run.delta,'FontWeight','bold','ForegroundColor','blue');
end
if (RunSet(k).RunStats.nIter >= str2double(get(gui.Term.nIter,'String')))
   set(Run.nIter, 'FontWeight','bold','ForegroundColor','red');
end
if (RunSet(k).RunStats.nFunc >= str2double(get(gui.Term.nFunc,'String')))
   set(Run.nFunc, 'FontWeight','bold','ForegroundColor','red');
end
if (RunSet(k).RunStats.time >= str2double(get(gui.Term.time,'String')))
   set(Run.time, 'FontWeight','bold','ForegroundColor','red');
end
if (RunSet(k).RunStats.nFails >= str2double(get(gui.Term.nFails,'String')))
   set(Run.nFails, 'FontWeight','bold','ForegroundColor','red');
end
if (RunSet(k).RunStats.stopRun)
   set(Run.stop, 'FontWeight','bold','ForegroundColor','red');
end
end

%*******************************************************************************
% solution_gui: View MADS optimal solution.
% ------------------------------------------------------------------------------
% Called by: results_gui (callback: View Solution pushbutton)
% Calls:     viewPrevious, viewNext
% VARIABLES:
%  iterate     = the best solution found by MADS
%    .x        =   vector of continuous variable values
%    .p        =   cell array of categorical variables
%  fig         = handle for the NOMADm GUI
%  gui.func    = handles for the NOMADm functions
%  gui_var     = handles for all GUI variables
%    .digits   =   number of digits accuracy to present on screen
%  gui.problem = name of current optimization problem
%  label       = text for "Best Feasible" or "Least Infeasible" solution
%  nX          = number of continuous variables
%  nP          = number of categorical variables
%  maxPage     = maximum number of variables displayed per page
%  nPages      = number of pages to display
%  onoff       = cell array containing "on" and "off"
%  leftmargin  = left margin position of continuous variables display
%  catmargin   = left margin position of categorical variables display
%  pos         = vertical position of next variable
%  xDisplay    = string conversion of the continuous variables
%  x(k)        = handles of the k-th displayed continuous variable
%    .ind      =   handle for the index label
%    .val      =   handle for the variable values
%  previousX   = handle for the continuous variables Previous push button
%  nextX       = handle for the continuous variables Next push button
%  pDisplay    = string conversion of the categorical variables
%  p(k)        = handles of the k-th displayed categorical variable
%    .ind      =   handle for the index label
%    .val      =   handle for the variable values
%  previousP   = handle for the categorical variables Previous push button
%  nextP       = handle for the categorical variables Next push button
%*******************************************************************************
function solution_gui(~,~,k,iterate) %#ok

fig      = findobj('Tag', 'NOMADm_GUI');
gui_var  = getappdata(fig,'gui_var');
gui      = getappdata(fig,'gui');
label    = get(gcbo,'UserData');
nX       = length(iterate.x);
nP       = length(iterate.p);
maxPage  = 14;
onoff    = {'on','off'};

% Set up "Display Solution" figure window
fig = figure('Name', ['MADS ',label,' Solution for Problem: ', ...
                      get(gui.problem, 'String'),', Run #' int2str(k)], ...
   'DefaultUIControlUnits',          'normalized', ...
   'DefaultUIControlFontUnits',      'normalized', ...
   'DefaultUIControlFontName',       'Helvetica',  ...
   'DefaultUIControlFontSize',        0.375, ...,
   'DefaultUIControlStyle',           'text', ...
   'Units',                           'normalized', ...
   'Position',                        [0.15 0.10 0.8 0.7], ...
   'MenuBar',                         'none', ...
   'NumberTitle',                     'off');

% Set display margins
if nP
   leftmargin = .08;
   catmargin  = .60;
else
   leftmargin = .30;
   catmargin  = .60;
end

% Display continuous variables Solution Labels
uipanel;
uicontrol(fig, 'String',             'CONTINUOUS VARIABLES', ...
               'FontWeight',         'bold', ...
               'HorizontalAlignment','center', ...
               'Position',           [leftmargin,.86,.34,.10]);
uicontrol(fig, 'String',             'Index', ...
               'HorizontalAlignment','left', ...
               'Position',           [leftmargin,.76,.07,.10]);
uicontrol(fig, 'String',             'Value', ...
               'HorizontalAlignment','right', ...
               'Position',           [leftmargin+.09,.76,.23,.10]);

% Display continuous variables
xDisplay = cell(maxPage,1);
for k = 1:nX
   xDisplay{k} = num2str(iterate.x(k),gui_var.Options.digits);
end

x(maxPage).ind = [];
x(maxPage).val = [];
for k = 1:maxPage
   pos = .76 - .04*k;
   x(k).ind = uicontrol(fig, ...
                 'String',             int2str(k), ...
                 'HorizontalAlignment','left', ...
                 'Position',           [leftmargin, pos, .05, .08]);
   x(k).val = uicontrol(fig, ...
                 'String',             xDisplay{k}, ...
                 'HorizontalAlignment','right', ...
                 'Position',           [leftmargin+.07, pos, .25, .08]);
end
set([x(nX+1:end).ind, x(nX+1:end).val]', 'Visible','off');

% Display Previous/Next pushbuttons for continuous variables
previousX = uicontrol(fig, ...
                      'Style',    'pushbutton', ...
                      'String',   'Previous',   ...
                      'Visible',  'off', ...
                      'Position', [leftmargin+.07, .05, .08, .06], ...
                      'Callback', {gui.func.changeView,gui_var.Options.digits});
nextX     = uicontrol(fig, ...
                      'Style',    'pushbutton', ...
                      'String',   'Next', ...
                      'Visible',  onoff{2 - (nX > maxPage)}, ...
                      'Position', [leftmargin+.19, .05, .08, .06], ...
                      'Callback', {gui.func.changeView,gui_var.Options.digits});
set(previousX,'UserData',{-1,x,iterate.x,nextX});
set(nextX,    'UserData',{ 1,x,iterate.x,previousX});

% If MVP problem:
if nP

   % Convert the categorical variables into a cell array of strings
   [pDisplay{1:maxPage}] = deal('');
   for k = 1:nP
      if ischar(iterate.p{k})
         pDisplay{k} = iterate.p{k};
      else
         pDisplay{k} = num2str(iterate.p{k});
      end
   end

   % Display categorical variables Solution Labels
   uicontrol(fig, 'String',             'CATEGORICAL VARIABLES', ...
                  'FontWeight',         'bold', ...
                  'HorizontalAlignment','center', ...
                  'Position',           [catmargin .86 .34 .10]);
   uicontrol(fig, 'String',             'Index', ...
                  'HorizontalAlignment','left', ...
                  'Position',           [catmargin .76 .07 .10]);
   uicontrol(fig, 'String',             'Value', ...
                  'HorizontalAlignment','right', ...
                  'Position',           [catmargin+.09 .76 .23 .10]);

   % Display categorical variables
   p(maxPage).ind = [];
   p(maxPage).val = [];
   for k = 1:maxPage
      pos = .76 - .04*k;
      p(k).ind = uicontrol(fig, ...
                  'String',             int2str(k), ...
                  'HorizontalAlignment','left', ...
                  'Position',           [catmargin pos .05 .08]);
      p(k).val = uicontrol(fig, ...
                  'String',             pDisplay{k}, ...
                  'HorizontalAlignment','right', ...
                  'Position',           [catmargin+.07 pos .25 .08]);
   end
   set([p(nP+1:end).ind, p(nP+1:end).val]', 'Visible','off');

   % Display Previous/Next pushbuttons for categorical variables
   previousP = uicontrol(fig, ...
                 'Style',    'pushbutton', ...
                 'String',   'Previous',   ...
                 'Visible',  'off', ...
                 'Position', [catmargin+.07, .05, .08, .06], ...
                 'Callback', {gui.func.changeView,gui_var.Options.digits});
   nextP     = uicontrol(fig, ...
                 'Style',    'pushbutton', ...
                 'String',   'Next', ...
                 'Visible',  onoff{2 - (nP > maxPage)}, ...
                 'Position', [catmargin+.19, .05, .08, .06], ...
                 'Callback', {gui.func.changeView,gui_var.Options.digits});
   set(previousP,'UserData',{-1,p,iterate.p,nextP});
   set(nextP,    'UserData',{ 1,p,iterate.p,previousP});
end
end

%*******************************************************************************
% changeView: Change the view of a multi-page list of variable values.
% ------------------------------------------------------------------------------
% Called by: solutions_gui (callback: Previous/Next pushbuttons)
% VARIABLES:
%   digits      = number of significant digits to display
%   UserData    = storage for key items needed
%   code        = an integer: magnitude = page #, sign = which button pressed
%   solution    = handles for the on-screen display of the variable values
%     .ind      =   handles for the variable indices
%     .val      =   handles for the variable values
%   iterate     = the solution vector
%   otherhandle = handle of the push button not pressed
%   n           = number of variables
%   whichButton = pushbutton indicator: Next = +1, Previous = -1
%   iPage       = display page number
%   maxPage     = maximum number of variables displayed on one page
%   lastPage    = the total number of display pages
%   newcode     = updated code variable assigned back to the pushbutton UserData
%   nDisplay    = number of variables displayed on the current page
%   ind         = index of variable to be displayed
%   value       = variable value converted to string for display
%*******************************************************************************
function changeView(~,~,digits) %#ok

% Retrieve pushbutton parameters
UserData = get(gcbo, 'UserData');
[code,solution,iterate,otherhandle] = deal(UserData{:});

% Set up page parameters and update pushbutton parameters for new page
n           = length(iterate);
whichButton = sign(code);
iPage       = abs(code);
iPage       = iPage + whichButton;
maxPage     = length(solution);
lastPage    = ceil(n/maxPage);
newcode     = iPage*whichButton;
set(gcbo,        'UserData',{ newcode,solution,iterate,otherhandle});
set(otherhandle, 'UserData',{-newcode,solution,iterate,gcbo});

% Set visibility of variable and index displays
set([gcbo; otherhandle], 'Visible', 'on');
nDisplay = maxPage;
if (iPage == 1)
   set(gcbo, 'Visible','off');
end
if (iPage == lastPage)
   nDisplay = mod(n,maxPage);
   if ~nDisplay, nDisplay = maxPage; end
   set(gcbo, 'Visible','off');
end

% Update display of indices and variables for new page
for k = 1:nDisplay
   ind = k + (iPage-1)*maxPage;
   if iscell(iterate)
      if ischar(iterate{ind})
         value = iterate{ind};
      else
         value = num2str(iterate{ind},digits);
      end
   else
      value = num2str(iterate(ind),digits);
   end
   set(solution(k).ind, 'Visible','on','String',int2str(ind));
   set(solution(k).val, 'Visible','on','String',value);
end
for k = nDisplay+1:maxPage
   set([solution(k).ind; solution(k).val], 'Visible','off');
end
end
