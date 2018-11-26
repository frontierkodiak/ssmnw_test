function [ ] = clusteranalysis_collate_Jacobians(timestamp, numReps, doDelete)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% ------------------ SETUP ------------------------------------------------

p = gcp

dataDir = ['/proj/ecoservice/users/x_alvcu/SSMNW/Results/', timestamp, '/'];
saveDir = ['/proj/ecoservice/users/x_alvcu/SSMNW/Results/', timestamp, '_collatedJacobians/'];

% ------------------ Loop over simulation result folders ------------------


% Checking necessary files are present.
psFile = dir([dataDir,'paramStruct_SSMNW.mat']);
jFiles = dir([dataDir,'Jacobian_*']);
if ( isempty(jFiles) && isempty(psFile) )
    fprintf('\n>>> NOTICE! >>>\n\n >>> Neither paramStruct nor Jacobians. Skipping dir: %s. >>>\n', dataDir)
    return
elseif ( isempty(jFiles) && ~isempty(psFile) )
    fprintf('\n>>> NOTICE! >>>\n\n >>> No Jacobians. Skipping dir: %s. >>>\n', dataDir)
    return
elseif ( ~isempty(jFiles) && isempty(psFile) )
    fprintf('\n>>> NOTICE! >>>\n\n >>> No paramStruct. Skipping dir: %s. >>>\n', dataDir)
    return
end

% Load parameter matrix
foo = load([dataDir,'paramStruct_SSMNW.mat']);
params = foo.params;
clear foo

% Define number of scenarios and replicates
numScens = size(params.data,1);
% numReps  = size(dir([dataDir,'Jacobian_scenID_', int2str(numScens), '*']),1);

% Bug check
if ( numScens*numReps ~= size(jFiles,1) )
    fprintf('\n>>> ERROR! >>>\n\n >>> File number incorrectly defined. >>>\n')
    doDelete = 0;
end


% ------------------ Loop over scenarios -----------------------------

tic
fprintf('\n>>> Starting loop for Folder: %s >>>\n>>> %i scenarios and %i replicates. >>>\n>>> Time: %s >>>\n', dataDir, numScens, numReps, datestr(now))

parfor sc = 1:numScens
    
    % Create matrix to hold Jacobians
    S = params.data(sc, params.S);
    repMat = zeros(S,S,numReps);
    
    % Loop over replicates and put individual Jacobians into repMat
    for r = 1:numReps
        fileName  = [dataDir,'Jacobian_scenID_', int2str(sc), '_replicate_', int2str(r), '.csv'];
        filecheck = dir(fileName);
        if ~isempty(filecheck)
            J = csvread(fileName);
            repMat(:,:,r) = J;
        else
            repMat(:,:,r) = NaN;
        end
    end %of replicate loop r=1:numReps
    
    % Save repMat, i.e. the collated Jacobians
    saveName = [saveDir, 'allJacobians_scenID_', num2str(sc), '.mat'];
    obj_parsave(saveName, repMat);
    
end %of scenario loop sc=1:numScens

% Delete individual Jacobian files (AFTER saving is safely done!)
if doDelete
    deleteFiles = dir([dataDir,'Jacobian_scenID_*']);
    if ~isempty(deleteFiles)
        cd(dataDir)
        delete(deleteFiles.name);
    end
end

fprintf('\n>>> Done with Folder: %s >>>\n>>> At: %s >>>\n', dataDir, datestr(now))
toc

end

