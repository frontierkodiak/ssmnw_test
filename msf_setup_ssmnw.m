nary function [ runProfiling, scurr, project, replicates, ...
    S, numGuilds, guildStruct, JacobianSum, intra, ...
    wGuildC, wGuildCompType, wGuildCompPairwise, ...
    wGuildCompSymmetry, wGuildCompCorr, wGuildWeightParams, ...
    btwGuildC, binNestLevel, quantNestLevel, numIntTypes,...
    intTypes, intTypeSymmetry, intTypeCorr, intWeightParams, ...
    validity, resMat, resMatFileName, ...
    repMat, JacobianFileName, saveJac, workspaceFileName ...
    ] = msf_setup_ssmnw( row, saveDir )
%
%MSF_SETUP_SSMNW Gives all input needed for main_ssmnw and its
%subfunctions. Here is where you set defaults, options, create data
%pathways and so on.
%
%   Included in the setup function is:
%
%   1) Setting options for running and saving of profiling (do or do not)
%   2) How to deal with warnings
%   3) Data paths to parameter matrix and where to save results.
%   4) Settings for random number generation.
%   5) Default values of scenario and model parameters.
%   6) The possibility of assigning values of model parameters using a
%       parameter matrix. If row = 0, then the default values are used. If
%       row > 0, the values on row "row" in the parameter matrix will be
%       used, for the parameters that are defined in the parameter matrix.
%       Parameters not defined in the parameter matrix keep their default
%       values.
%   7) Scenario validity checks
%   8) Create structures for saving.
%   9) Final operations on parameters.
%
%   @INPUT
%   row     Scenario, i.e. row index in parameter matrix. If row=0, the
%           default values hard coded in this file will be used.
%
%   @OUTPUT
%   runProfiling        Logical scalar. Whether to run profiling or not.
%   saveDir             String. Directory name, where to save results and all.
%   scurr               Struct. Random generator seed.
%   project             String. Name of sub project.
%   trials              Integer. Number of trials or replicates
%   S                   Integer. Species richness (total)
%   numGuilds           Integer or numerical vector. Number of Guilds                            [1, 2, 3]
%   guildStruct         Presence of Guild structure                 [0, 1], i.e. 0=No, 1=Yes
%   JacobianSum         Scenario independent Jacobian sum or not    [0, 1], i.e. 0=Jacobian sum is not held constant, 1=Jacobian sum is held constant across scenarios with the same S and C
%   intra               Intraspecific competition.                  Integer.
%   wGuildC             Connectance of Within-Guild interactions    Fraction btw 0 and 1.
%   wGuildCompType      Type of Within-Guild competition            [1, 2], i.e. 1='diffuse', 2='diffuseRand'
%   wGuildCompPairwise  Whether binary wG ints are to be pairwise   [0, 1], i.e. 0=No (aij=1 is independent of aji=1), 1=Yes (if aij=1, aji=1 too)
%   wGuildCompSymmetry  Symmetry of wg competition                  [0, 1, 2, 3], i.e. 0=independence, 1=perfect symmetry, 2=specified numerical correlation, 3=rank-order correlation. 
%   wGuildCompCorr     	Correlation of pairwise wg IS               Value btw -1 and 1. The desired Pearson's correlation coeff. Parameter used when intTypeSymmetry == 2.
%   wGuildWeightParams  Within-Guild IS distribution params         [Mean, STD] of a normal distribution
%   btwGuildC           Connectance of Btw-Guild interactions       Fraction btw 0 and 1.
%   binNestLevel        Level of binary nestedness                  Value btw -1 and 1 (-1 = perfect antinestedness, 0 = random, 1 =  perfect nestedness, other values are in between)
%   quantNestLevel      Level of quantitative nestedness            Value btw -1 and 1 (-1 = perfect antinestedness, 0 = random, 1 =  perfect nestedness, other values are in between)
%   numIntTypes         Number of Btw-Guild Interaction Types       [1, 2, 3]
%   intTypes            The names of the Btw-Guild Int Types        [1, 2, 3], i.e. 1='competition', 2='mutualism', 3='trophic'
%   intTypeSymmetry     Symmetry of pairwise Btw-Guild interactions [0, 1, 2, 3], i.e. 0=independence, 1=perfect symmetry, 2=specified numerical correlation, 3=rank-order correlation. 
%   intTypeCorr         Correlation of pairwise Btw-Guild IS.       Value btw -1 and 1. The desired Pearson's correlation coeff. Parameter used when intTypeSymmetry == 2.
%   intWeightParams     Between-Guild int strength distribution     [Mean, STD] of a normal distribution
%
%   OTHER NOTABLE VARIABLES
%   wgcIgnore           Ignore wGuildC; instead use btwGuildC       [0, 1], i.e. 0=No, 1=Yes
%
%   @AUTHORS
%   Alva Curtsdotter, Post doc, Dep of Ecology, SLU Ultuna, 2015-04-14

% 1) RUN PROFILING -----------------------------------------------------------

runProfiling = 0;                                                           % Set to 0 if you do not wish to run the profiler!

if runProfiling
    profile on
end


% 2) TURN OFF WARNINGS -------------------------------------------------------

warning('off','MATLAB:HandleGraphics:noJVM')                                % Turn off warning generated from figure() in an -nojvm environment.


% 3) DATA PATHS & SAVING OPTIONS ---------------------------------------------

% Parameter matrix directory paths
if isunix                           % To run on Linux platform, i.e. the cluster.
    PS_DIR = '/scratch/local/';
elseif ispc                         % To run on Windows platform, i.e. my computer   
    addpath('C:\Users\icurtsd\Documents\Project Atlanta\Project Structure and stability of mutualistic networks\code\')
    PS_DIR = 'C:\Users\icurtsd\Documents\Project Atlanta\Project Structure and stability of mutualistic networks\code\';
elseif ismac
    fprintf('\n>>> NOTICE! >>>\n\n >>> In setup_ssmnw:\n >>> Mac? Must be Connor and Loy.\n Hey guys, make sure the all directory paths are written in a mac-friendly way, ok?\n It might differ from PC, but I dont know. >>>\n')
else
    fprintf('\n>>> WARNING! >>>\n\n >>> In setup_ssmnw:\n >>> Unknown platform. No PS_DIR defined. >>>\n')
end

% % Save directory path
% saveDir = 'C:\Users\icurtsd\Documents\Project Atlanta\Project Structure and stability of mutualistic networks\data\main_test\';
% timeStamp = datestr(now,'yyyymmmddTHHMM');
% saveDir = [saveDir, 'testresults_', timeStamp, '\'];
% if 7 ~= exist(saveDir, 'dir')
%     mkdir(saveDir)
% end


% 4) SETTINGS FOR RANDOM NUMBER GENERATION -----------------------------------

rng('shuffle')
scurr = rng;


% 5) MODEL CONSTANTS ---------------------------------------------------------

% Technical parameters
project = 1;                            % Name of sub project. [1, 2], i.e. 1='guild', 2='nestedness'
trials  = 10;                           % Number of trials or replicates; PC setup. Integer.
saveJac = trials;                       % Number of Jacobians to save. Integer.
                
% Whole network  parameters
S           = 20;                       % Species richness (total). Even integer.
numGuilds   = 2;                        % Number of Guilds. [1, 2, 3]
guildStruct = 1;                        % Presence of Guild structure. [0, 1], i.e. 0=No, 1=Yes
JacobianSum = 0;                        % Scenario independent Jacobian sum or not. [0, 1], i.e. 0=Jacobian sum is not held constant, 1=Jacobian sum is held constant across scenarios with the same S and C

% Between guild interaction structure
btwGuildC       = 0.5;                  % Connectance of Btw-Guild interactions. Fraction btw 0 and 1.
binNestLevel    = 1;                    % Binary nestedness. 1=perfectly nested, 0=random, -1="perfectly" antinested
quantNestLevel  = 1;                    % Quantitative nestedness. 1=perfectly nested, 0=random, -1=perfectly antinested (sensu Staniczenko, Ghost paper)
numIntTypes     = 1;                    % Number of Btw-Guild Interaction Types. [1, 2, 3]
intTypes        = 3;                    % The names of the Btw-Guild Int Types. [1, 2, 3], i.e. 1='competition', 2='mutualism', 3='trophic'
intTypeSymmetry = 3;                    % Symmetry of pairwise Btw-Guild interactions. [0, 1, 2, 3], i.e. 0=independence, 1=perfect symmetry, 2=specified numerical correlation, 3=rank-order correlation. 
intTypeCorr     = 0.7;                  % Correlation of pairwise Btw-Guild interaction strengths when intTypeSymmetry=2.
intWeightMean   = 0;                    % Btw-Guild IS distribution params. Mean of a normal distribution
intWeightStd    = sqrt(1/4);            % Btw-Guild IS distribution params. STD of a normal distribution

% Within guild interaction structure
intra              = -1;                % Intraspecific competition. Negative integer.
wGuildC            = 0;                 % Connectance wGuildL/(S/2)^2 of Within-Guild interactions. Fraction btw (S/2)/(S/2)^2 (i.e. intra) and 1.
wgcIgnore          = 0;                 % If 1, set wGuildC to btwGuildC
wgcSubstitutive    = 1;                 % If 0, additive design. If 1, substitutive design.
wGuildCompType     = 2;                 % Type of Within-Guild competition. [1, 2], i.e. 1='diffuse', 2=other (formerly 'diffuseRand', which now corresponds to wGuildCompType=2 and wGuildCompSymmetry=0)
wGuildCompPairwise = 1;                 % Whether the within-guild interactions are pairwise (1) or not (0).
wGuildCompSymmetry = 0;                 % Symmetry of pairwise Within-Guild competition. [0, 1, 2, 3], i.e. 0=independence, 1=perfect symmetry, 2=specified numerical correlation, 3=rank-order correlation. 
wGuildCompCorr     = 0.75;              % Correlation of pairwise Within-Guild competition strengths when wGuildCompSymmetry=2.
wGuildWeightMean   = 0;                 % Within-Guild IS distribution params. Mean of a normal distribution
wGuildWeightStd    = sqrt(1/4);%sqrt(1/16);% Within-Guild IS distribution params. STD of a normal distribution
wgWeightIgnore     = 1;                 % If 1, set wGuildWeights to btwGuildWeights (i.e. intWeight)

% 6) MODEL VARIABLES ---------------------------------------------------------
% Any of the above constants can be made into scenario variables, if
% 1) they are part of the parameter matrix, and
% 2) the code is written to accomodate variation in the parameter in question.
% The parameters that are present in the parameter matrix will have their
% values above overwritten by the values in the parameter matrix.

if row == 0
    fprintf('\n>>> NOTICE! >>>\n\n >>> In setup_ssmnw:\n >>> No parameter matrix. Hard coded default values are used. >>>\n')
else
    load([PS_DIR 'paramStruct_SSMNW.mat'])  % Structure called params is loaded.
    for i = 1:numel(params.names)
        column = eval(['params.',params.names{i}]);
        eval([sprintf(params.names{i}), '=', num2str(params.data(row, column))]);
    end
end                  

% 7) CHECK SCENARIO VALIDITITY -----------------------------------------------

nSinG     = round(S/2);                        % Number of species in each guild.
nlinksBtwG = round(btwGuildC*nSinG^2);         % Number of links between guilds.

% Check that the minimum requirement of 1 btw guild link per species is met.
if nlinksBtwG < nSinG
    validity = 0;
else
    validity = 1;
end


% 8) PREPARING SAVING STRUCTURES ---------------------------------------------

if validity == 1
    resMat = zeros(trials,3);    % Structure for saving results for local stability, resilience, and reactivity.
    repMat = zeros(S,S,saveJac); % Structure for saving the Jacobians.
else
    resMat = NaN(trials,3);
    repMat = [];
end

if row == 0
    scenarioName = '_scenID_default';
else
    scenarioName = ['_scenID_', int2str(row)];    
end

paramsFileName     = [ saveDir, 'paramStruct_SSMNW.mat'];
resMatFileName     = [ saveDir, 'ResMat', scenarioName, '.mat'];
JacobianFileName   = [ saveDir, 'allJacobians', scenarioName, '.mat'];
workspaceFileName  = [ saveDir, 'WorkSpaceParams', scenarioName, '.mat'];

% Saving the parameter matrix to the saveDir.
if ( row == 1 )
    save(paramsFileName, 'params')
end
    
% 9) SOME FINAL OPERATIONS ON VARIABLES TO SHAPE THEM OR CREATE NEW ONES -----

% Network parameters
wGuildWeightParams = [wGuildWeightMean, wGuildWeightStd];
intWeightParams    = [intWeightMean, intWeightStd];

if wgWeightIgnore
   wGuildWeightParams = intWeightParams;
end

if wgcIgnore
   wGuildC = btwGuildC; 
end

if ( wgcSubstitutive && wgcIgnore )
    btwGuildC = btwGuildC/2;    
    wGuildC   = btwGuildC;
end

% Technical parameters
replicates = 1:trials;   


end

