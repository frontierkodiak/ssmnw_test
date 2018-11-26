function [ ] = analysis_collate_Jacobians()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% ------------------ SETUP ------------------------------------------------

% Define simulation result folders (dirs) to go through
mainDir = 'C:\Users\icurtsd\Documents\Project Atlanta\Project Structure and stability of mutualistic networks\data\main_test\';
dirs = dir([mainDir, 'testresults_*']);


% ------------------ Loop over simulation result folders ------------------

for d = 1:length(dirs)
    
    % Define directory
    dataDir = [mainDir, dirs(d).name, '\'];
%     saveDir = dataDir;
    
    % Checking necessary files are present.
    psFile = dir([dataDir,'paramStruct_SSMNW*']);
    jFiles = dir([dataDir,'Jacobian_*']);
    if ( isempty(jFiles) && isempty(psFile) )
        fprintf('\n>>> NOTICE! >>>\n\n >>> Neither paramStruct nor Jacobians. Skipping dir: %s. >>>\n', dataDir)
        continue
    elseif ( isempty(jFiles) && ~isempty(psFile) )
        fprintf('\n>>> NOTICE! >>>\n\n >>> No Jacobians. Skipping dir: %s. >>>\n', dataDir)
        continue
    elseif ( ~isempty(jFiles) && isempty(psFile) )
        fprintf('\n>>> NOTICE! >>>\n\n >>> No paramStruct. Skipping dir: %s. >>>\n', dataDir)
        continue
    end
    
    % Load parameter matrix
    load([dataDir,psFile.name]);
 
    
    % Define number of scenarios and replicates
    numScens = size(params.data,1);
    numReps  = size(dir([dataDir,'Jacobian_scenID_', int2str(numScens), '*']),1);
    
    % Bug check
    if ( numScens*numReps ~= size(jFiles,1) )
        fprintf('\n>>> ERROR! >>>\n\n >>> File number incorrectly defined. >>>\n')
        % Give keyboard control if deviations occur.
        keyboard
        % Go to next turn in folder loop if deviations occur.
        % continue
    end    
    
    
    % ------------------ Loop over scenarios -----------------------------
        
    tic
    fprintf('\n>>> Starting loop for Folder: %s >>>\n>>> %i scenarios and %i replicates. >>>\n>>> Time: %s >>>\n', dirs(d).name, numScens, numReps, datestr(now))
    
    missing = zeros(numScens*numReps - size(jFiles,1),2);
m = 1;

for sc = 1:numScens
    
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
            missing(m,1)  = sc;
            missing(m,2)  = r;
            m = m + 1;
        end
    end %of replicate loop r=1:numReps
    
    % Save repMat, i.e. the collated Jacobians
    saveName = [dataDir, 'allJacobians_scenID_', num2str(sc), '.mat'];
    save(saveName, 'repMat')
    
end %of scenario loop sc=1:numScens
  
    % Delete individual Jacobian files (AFTER saving is safely done!)
%     keyboard
    deleteFiles = dir([dataDir,'Jacobian_scenID_*']);
    if ~isempty(deleteFiles)
        cd(dataDir)
        delete(deleteFiles.name)
        clear jFiles
    end
    
    fprintf('\n>>> Done with Folder: %s >>>\n>>> At: %s >>>\n', dirs(d).name, datestr(now))
    toc
    
    
end %of dataDir loop

end

