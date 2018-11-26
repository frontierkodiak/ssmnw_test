function [ ] = wrapper_create_parameter_matrix()
%
% WRAPPER FOR CREATE_PARAMETER_MATRIX. Specify parameters and their levels
% and create the parameter matrix.
%
% This wrapper runs function create_parameter_matrix. A fully factorial 
% numerical matrix will be created, where each row corresponds to a
% scenario to be run and each row thus has a unique combination of scenario
% parameter values. This matrix and the paramater names (and the column
% number on which this variable occurs in the matrix) are stored in the
% struct params. You define the parameter names and the levels which they'r
% supposed to take in this current function. If you want combinations of
% parameter levels that are not fully factorial you will need to prune the
% matrix (section 8). You will need to specify project, name the variables,
% explore and/or define the variable levels, provide them in the function
% call, specify the save directory, and create/apply pruning options
% manually in the appropriate section. 
%
% Function structure
% 1) Define project
% 2) Names of the variables
% 3) Explore potential variable levels
% 4) Define variable levels
% 5) Make sanity check of chosen levels (disabled)
% 6) Where to save the resulting parameter matrix.
% 7) Run create_parameter_matrix
% 8) How to prune the parameter matrix
%
%   @AUTHORS
%   Alva Curtsdotter, Postdoc, Emory University, 2018-05-31
%   alva.curtsdotter@emory.edu
%
% -------------------------------------------------------------------------

% List of (potential) scenario variables, with explanation
%
% Variable          Explanation                                 Possible values
% project           Name of sub project.                        [1, 2], i.e. 1='guild', 2='nestedness'
% S                 Species richness (total)                    Integer
% numGuilds         Number of Guilds                            [1, 2, 3]
% guildStruct       Presence of Guild structure                 [0, 1], i.e. 0=No, 1=Yes
% wGuildC           Connectance of Within-Guild interactions    Fraction btw 0 and 1.
% wGuildCompType    Type of Within-Guild competition            [1, 2], i.e. 1='diffuse', 2='diffuseRand'
% wGuildWeightStd   Standard dev. of WG interaction strength    Integer.
% btwGuildC         Connectance of Btw-Guild interactions       Fraction btw 0 and 1.
% numIntTypes       Number of Btw-Guild Interaction Types       [1, 2, 3] 
% intTypes          The names of the Btw-Guild Int Types        [1, 2, 3], i.e. 1='competition', 2='mutualism', 3='trophic'
% intTypeSymmetry   Symmetry of pairwise Btw-Guild interactions [0, 1], i.e. 0=independence, 1=symmetry 
% binNestLevel      Level of binary nestedness                  Numerical btw -1 and 1. -1: antinested, 0: intermediate/random, 1: perfectly nested
% quantNestLevel    Level of quantitative nestedness            Numerical btw -1 and 1. -1: antinested, 0: intermediate/random, 1: perfectly nested
% JacobianSum       Scenario independent Jacobian sum or not    [0, 1], i.e. 0=Jacobian sum is not held constant, 1=Jacobian sum is held constant across scenarios with the same S and C
% trials            Number of trials or replicates              Integer
% saveJac           Number of Jacobians to save.                Integer.
%
% -------------------------------------------------------------------------


% 1) Define project! ------------------------------------------------------

project         =  1;                         % 1='guild', 2='nestedness'


% 2) Names of the variables or building the parameter matrix ---------------

%   Make sure the names of the variables are exactly the same as in the
%   function in which the parameter matrix is to be used.   

if project == 1
    variableNames = {'project', 'S', 'guildStruct',...
        'wGuildC', 'wgcIgnore', 'wgcSubstitutive', ...
        'btwGuildC', 'intTypes', 'intTypeSymmetry',  'intWeightStd', ...
        'trials', 'saveJac'};
elseif project == 2
    variableNames = {'project', 'S', ...
        'wGuildC', 'wgcIgnore', 'wgcSubstitutive',  ...
        'btwGuildC', 'intTypes', 'intTypeSymmetry', 'intWeightStd', ...
        'binNestLevel', 'quantNestLevel', ...
        'trials', 'saveJac'};
end


% 3) Explore variable levels -----------------------------------------------

% Play around a bit if u want, to help select levels of S and between guild C.
if project == 1
    S               = [250];
    s = S/2;
    minCp1a = 1./s; % For reference, will give you minimum C for your S, for project 1, when additive design.
    minCp1s = 2./s; % For reference, will give you minimum (btwGuild)C for your S, for project, when substitutive design (i.e minC for wGuildC=0, so that when you include wGuildC>0
    minCp2 = (S+2)./(s).^2; % For reference, will give you minimum C for your S, for project 2.
    for i = 1:length(S)
       maxCp2(i) = ( 2 - sum(0:(s(i)-3)) + ( s(i)-1 )^2 ) / s(i)^2; 
       maxCp1(i) = s(i)*(s(i)-1)/s(i)^2;
    end
    foo=[S; minCp1a; minCp1s; minCp2; maxCp2; maxCp1]
    
    btwGClevelsSub = 5;
%     btwGClevelsAdd = 3;
    btwGCVecSub = zeros(btwGClevelsSub,length(S));
%     btwGCVecAdd = zeros(btwGClevelsAdd,length(S));
    for whatever = 1:length(S)
        % Old something
        %         btwGCVec(:,whatever) = foo(3,whatever):(foo(4,whatever)-foo(3,whatever))/(btwGClevels-1):foo(4,whatever);
        
        % Substitutive design
        btwGCVecSub(:,whatever) = foo(3,whatever):(foo(6,whatever)-foo(3,whatever))/(btwGClevelsSub-1):foo(6,whatever);
        
%         % Additive design
%         btwGCVecAdd(:,whatever) = foo(2,whatever):(foo(6,whatever)-foo(2,whatever))/(btwGClevelsAdd-1):foo(6,whatever);
    end
%     btwGClevels = 3+1;
%     btwGCVec = [minCp1; btwGCVec];
% btwGCVec = btwGCVecSub;
btwGCVec = [ minCp1a; btwGCVecSub];
    allBtwGuildCLevels = unique((sort(btwGCVec(:))))';

elseif project == 2
    S               = 100; %[10 16 32 64 128];%[10 16 32 64 128];%64;%[10 20 128];
    s = S/2;
    minCa = (S+2)./(s).^2;
    minCs = 2*(S+2)./(s).^2;
    clear maxC
    for i = 1:length(S)
       maxC(i) = ( 2 - sum(0:(s(i)-3)) + ( s(i)-1 )^2 ) / s(i)^2; 
    end
    foo=[S; minCa; minCs; maxC]
    
    btwGClevels = 3;
    btwGCVec = zeros(btwGClevels,length(S));
    for whatever = 1:length(S)
        btwGCVec(:,whatever) = foo(3,whatever):(foo(4,whatever)-foo(3,whatever))/(btwGClevelsSub-1):foo(4,whatever);
    end
    links = round(btwGCVec.*(repmat(foo(1,:),btwGClevels, 1)/2).^2)     
    allBtwGuildCLevels = (sort(btwGCVec(:)))';
end

% Play around a bit if u want, to help select levels of interaction
% strengths.
if project == 1
    jacobianC       = 0.25; %                          % Connectance of Jacobian (L/S^2, but counting only off-diag elements in L) 
    
    intTypes        = [2, 3];%3;%
    guildStruct     = [0, 1];
    wgcIgnore       = [0, 1];
    fooNames        = {'guildStruct', 'wgcIgnore', 'intTypes'};    
    tmpParams = create_parameter_matrix([], '', fooNames, ...
       guildStruct, wgcIgnore, intTypes);
   
%     x = [0.15, 0.3; 0.15, 0.3; 0.45, 0.75; 0.45, 0.75; 3, 15; 1, 15; 1, 2.75; 1, 2.75];    
%     x = [0.1, 15;0.1, 15;0.1, 15;0.1, 15];
%     x = [0.15, 0.3; 0.15, 0.3; 0.45, 0.75; 0.45, 0.75; 0.1, 12; 0.1, 12; 0.1, 2.75; 0.1, 2.75];     
%     x = [0.1, 0.2; 0.1, 0.2; 0.4, 0.8; 0.4, 0.8; 2.2, 3.2; 1.8, 2.8; 0.8, 1.3; 0.8, 1.3  ]; +
%     x = [ 0.1, 12; 0.1, 12; 0.1, 12; 0.1, 12; 0.1, 12; 0.1, 12; 0.1, 12; 0.1, 12 ];
     x = [ 0.01, 0.25; 0.01, 0.25; 0.2, 1; 0.2, 1; 2, 12; 2, 12; 0.95, 2; 0.95, 2 ];
    xRange = linspace(0,1,50);
    xRange = x(:,1) + xRange.*(x(:,2) - x(:,1));
    fooNames
    [tmpParams.data, x]
    
    sigmaRange = xRange./repmat(sqrt(S*jacobianC),size(xRange,1),size(xRange,2));    
    allSigmas = unique(sigmaRange)';
    sigmaTable = [tmpParams.data, sigmaRange];
end


% 4) Define the levels. Note that S and C can be set above.  ------

if project == 1
%     S               = [8 32 128];
%     numGuilds       = 2;
    guildStruct     = [0, 1];                     % 0=No, 1=Yes
    wGuildC         = 0; %[0 allBtwGuildCLevels];   %[0 0.5 1];                  % minimum and maximum C values that we are interested in. Must check which are possible given scenario later...
    wgcIgnore       = [0, 1];
    wgcSubstitutive = 1;
%     wGuildCompType  = 2;                            %  1='diffuse', 2='diffuseRand'
%     wGuildWeightStd = [sqrt(1/16), sqrt(1/4)];
    
    btwGuildC       = 2*jacobianC;%true for substitutive design w wgC=0; when wgC=bgC, bgC=jacC . %allBtwGuildCLevels;% [0.1 0.25 0.5 0.95];                  % minimum and maximum C values that we are interested in. Must check which are possible given scenario later...
%     numIntTypes     = 1;
    intTypes        = [2, 3];%3;%  2;% [2, 3];% [1, 2, 3]; % 1='competition', 2='mutualism', 3='trophic'
    intTypeSymmetry = [0, 0.25, 0.5, 0.75, 1]; % 0=independence, 1=symmetry
    intWeightStd    =  allSigmas;% [sqrt(1/64), sqrt(1/32), sqrt(1/16)];
%     JacobianSum     = 0; % 0=Jacobian sum is not held constant, 1=Jacobian sum is held constant across scenarios with the same S and C
    trials          = 100;
    saveJac         = 10; 
    
elseif project == 2
    
%     S               = [10 16 22 64 88 120];% [10 20 30 128]; %[8:2:128];
%     numGuilds       = 2;
%     guildStruct     = 1;                          % 0=No, 1=Yes
    wGuildC         =  0;% [0 1];                   % minimum and maximum C values that we are interested in. Must check which are possible given scenario later...
    wgcIgnore       = [0, 1];
%     wGuildCompType  = 2;                            %  1='diffuse', 2='diffuseRand'
    wGuildWeightStd = [sqrt(1/16), sqrt(1/4)];%
    btwGuildC       = 0.3; %will be halved with sub design. with add des bgC doesnt change but jacC doubles % allBtwGuildCLevels; %[0.1 0.25 0.5 0.95];                  % minimum and maximum C values that we are interested in. Must check which are possible given scenario later...
%     numIntTypes     = 1;
    intTypes        = [2, 3]; % 1='competition', 2='mutualism', 3='trophic'
    intTypeSymmetry = 0; % 1=perfect symmetry, 0=independence (of values drawwn, but note that paiwise IS will be rankorder correlated by default). 
    binNestLevel    = [1, 0, -1];
    quantNestLevel  = [1, 0, -1];
%     JacobianSum     = 0; % 0=Jacobian sum is not held constant, 1=Jacobian sum is held constant across scenarios with the same S and C
    trials          = 100;
    saveJac         = 10; 
end


% % 5) Check sanity of variable combos ------------------------------------
% 
% wGuildL = repmat(wGuildC, length(S), 1).*(repmat((S/2).^2, length(wGuildC), 1)');
% if any(any(wGuildL<(S/2)' == 1 & repmat(wGuildC, length(S), 1) > 0 ))
%     fprintf('\n>>> NOTE! \n>>> Less within guild interactions than species in one or more scenarios. Decide how to deal with it!')
%     keyboard       
% end
% btwGuildL = repmat(btwGuildC, length(S), 1).*(repmat((S/2).^2, length(btwGuildC), 1)');
% if any(any(btwGuildL<(S/2)' == 1 & repmat(btwGuildC, length(S), 1) > 0 ))
%     fprintf('\n>>> NOTE! \n>>> Less between guild interactions than species in one or more scenarios. Decide how to deal with it!')
%     keyboard   
% end 
% 


% 6) Define directory where to save the parameter matrix s------------------

SAVE_DIR = 'C:\Users\Caleb\Desktop\Research\Nestedness Matrices\';
saveName = 'paramStruct_SSMNW.mat';


% 7) Run CREATE_PARAMETER_MATRIX. -----------------------------------------

%   Make sure the vectors with the variable levels are put into the
%   function call in the EXACT same order as they are in the variableNames 
%   array. This is VERY IMPORTANT.

if project == 1
    
    params = create_parameter_matrix(SAVE_DIR, saveName, variableNames, ...
        project, S, guildStruct,...
        wGuildC, wgcIgnore, wgcSubstitutive, ...
        btwGuildC, intTypes, intTypeSymmetry, intWeightStd, ...
        trials, saveJac);
    
elseif project == 2
    
    params = create_parameter_matrix(SAVE_DIR, saveName, variableNames, ...
        project, S, numGuilds, guildStruct,...
        wGuildC, wgcIgnore, wGuildCompType, wGuildWeightStd, ...
        btwGuildC, numIntTypes, intTypes, intTypeSymmetry, binNestLevel, quantNestLevel, ...
        JacobianSum, trials, saveJac);
end


% 8) IF YOU NEED TO PRUNE YOUR PARAMETER MATRIX ---------------------------

% Reload the parameter matrix (or structure actually) that you just created
% load(saveName)
params

% Keep only the "right" S-specific btwGuildC levels for each S level

% for whatever = 1:length(S)    
%     params.data = params.data( ~(params.data(:,params.S) == S(whatever) & ~ismember(params.data(: ,params.btwGuildC), btwGCVec(:,whatever))) ,:);
% %     params.data = params.data( ~(params.data(:,params.S) == S(whatever) & ~ismember(params.data(: ,params.wGuildC), [0;btwGCVec(:,whatever)] )) ,:);
% end

% Prune species richness levels depedening on intTypes

% params.data = params.data( ~(params.data(:,params.intTypes) == 1 & params.data(: ,params.S)>20) ,:);
% params.data = params.data( ~(params.data(:,params.intTypes) == 2 & params.data(: ,params.S)>20) ,:);
% params.data = params.data( ~(params.data(:,params.intTypes) == 3 & ismember(params.data(: ,params.S), [10, 14])  ) ,:);

% Prune away multiple levels of withinGuildIS for cases with wGuildC == 0

% params.data = params.data( ~(params.data(:,params.wGuildC) == 0 & params.data(:,params.wgcIgnore) == 0 & ~ismember(params.data(: ,params.wGuildWeightStd), wGuildWeightStd(1))) ,:);

% Only need one level of wGuildC when wgcIgnore==1

% params.data = params.data( ~(params.data(:,params.wGuildC) ~= wGuildC(1)  & params.data(: ,params.wgcIgnore)==1) ,:);

for therow = tmpParams.data(:,1)' 
    trueRows = params.data(:,params.guildStruct) == tmpParams.data(therow,tmpParams.guildStruct) & params.data(:,params.wgcIgnore) == tmpParams.data(therow,tmpParams.wgcIgnore) &  params.data(:,params.intTypes) == tmpParams.data(therow,tmpParams.intTypes);
    params.data = params.data( ~( trueRows & ~ismember(params.data(: ,params.intWeightStd), sigmaTable(therow,5:end)) ) ,:);
%     params.data = params.data( ~(params.data(:,params.S) == S(whatever) & ~ismember(params.data(: ,params.wGuildC), [0;btwGCVec(:,whatever)] )) ,:);
end


% Do some sanity checks

% unique(params.data(params.data(:,params.intTypes)==1,params.S))
% unique(params.data(params.data(:,params.intTypes)==2,params.S))
% unique(params.data(params.data(:,params.intTypes)==3,params.S))

unique(params.data(params.data(:,params.wgcIgnore)==1,params.wGuildC))
unique(params.data(params.data(:,params.wgcIgnore)==1,params.btwGuildC))
unique(params.data(params.data(:,params.wgcIgnore)==1,params.intWeightStd))
% unique(params.data(params.data(:,params.wgcIgnore)==1,params.wGuildWeightStd))

unique(params.data(params.data(:,params.wgcIgnore)==0,params.wGuildC))
unique(params.data(params.data(:,params.wgcIgnore)==0,params.btwGuildC))
unique(params.data(params.data(:,params.wgcIgnore)==0,params.intWeightStd))
% unique(params.data(params.data(:,params.wgcIgnore)==0,params.wGuildWeightStd))
% unique(params.data(params.data(:,params.wgcIgnore)==0 & params.data(:,params.wGuildC)==1,params.wGuildWeightStd))
% unique(params.data(params.data(:,params.wgcIgnore)==0 & params.data(:,params.wGuildC)~=1,params.wGuildWeightStd))
% unique(params.data(params.data(:,params.wgcIgnore)==0 & params.data(:,params.wGuildC)==0,params.wGuildWeightStd))

for therow = tmpParams.data(:,1)' 
    trueRows = params.data(:,params.guildStruct) == tmpParams.data(therow,tmpParams.guildStruct) & params.data(:,params.wgcIgnore) == tmpParams.data(therow,tmpParams.wgcIgnore) &  params.data(:,params.intTypes) == tmpParams.data(therow,tmpParams.intTypes);
    disp(length(params.data( trueRows,params.intWeightStd)))
%     params.data = params.data( ~(params.data(:,params.S) == S(whatever) & ~ismember(params.data(: ,params.wGuildC), [0;btwGCVec(:,whatever)] )) ,:);
end

% Make sure scenario id numbers are continous
params.data(:,1) = 1:size(params.data,1);

% Check it
params

% Save the pruned parameter matrix struct.
save(saveName, 'params')

end




