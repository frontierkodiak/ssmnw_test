function [ ] = main_ssmnw( row, saveDir )
%
%   @INPUT
%   row       Scenario number = Row to be read from the parameter matrix.
%   saveDir   Where to save output
%
%   @OUTPUT (data saved to disk)
%
%   @AUTHORS
%   Alva Curtsdotter, Post doc @ BrosiLab, Dep of Environmental Sciences,
%   Emory University, Atlanta, Georgia, USA. Code initiated 2018-06-05.
%
% -------------------------------------------------------------------------

% Temp for dev purposes
% row=0; saveDir=[];

tic
fprintf('\n\n>>> RUNNING SCENARIO NR %i <<<\n', row)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SETUP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n>>> Running Setup <<<\n')

[ runProfiling, scurr, project, replicates, ...
    S, numGuilds, guildStruct, JacobianSum, intra, ...
    wGuildC, wGuildCompType, wGuildCompPairwise, ...
    wGuildCompSymmetry, wGuildCompCorr, wGuildWeightParams, ...
    btwGuildC, binaryNestednessLevel, quantitativeNestednessLevel, ...
    numIntTypes, intTypes, intTypeSymmetry, intTypeCorr, intWeightParams, ...
    valid, resMat, resMatFileName, ...
    repMat, JacobianFileName, saveJac, workspaceFileName ] = ...
    msf_setup_ssmnw( row, saveDir );                                % Set model parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM STARTS HERE   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if valid
    
    for t = replicates                                                      % Loop over all the replicates for the chosen scenario.
        % t=1; %temporary for development purposes
        
        fprintf('\n\n>>> RUNNING REPLICATE NR %i <<<\n', t)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%   COMMUNITY STRUCTURE AND PARAMETARIZATION   %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % CREATE BINARY BETWEEN-GUILD INTERACTION STRUCTURE -----------------------
        %  this inlcludes creating the weight order
        
        if project == 1
            
            [ APchunk, PAchunk, nlinksBtwG ] = ...
                msf_create_binary_btwGuildStructure_ssmnw( ...
                project, S, numGuilds, btwGuildC, intTypeSymmetry );
            
        elseif project == 2
            
            [ APchunk, PAchunk, nlinksBtwG ] = ...
                msf_create_binary_btwGuildStructure_ssmnw( ...
                project, S, numGuilds, btwGuildC, intTypeSymmetry, ...
                binaryNestednessLevel, quantitativeNestednessLevel );
            
        end
     
        
        % CREATE WEIGHTED BETWEEN-GUILD INTERACTION STRUCTURE ---------------------
        % Generate random weights and distribute them according to order
        
        [ APchunk, PAchunk ] = ...
            msf_create_weighted_btwGuildStructure_ssmnw( ...
            APchunk, PAchunk, intWeightParams, nlinksBtwG, numIntTypes, ...
            intTypeSymmetry, intTypeCorr, intTypes );
        

        % CREATE WEIGHTED WITHIN-GUILD INTERACTION STRUCTURE ----------------------
        % Generate and distribute random weights
        
        [ PPchunk, AAchunk ] = ...
            msf_create_weighted_wGuildStructure_ssmnw( ...
            APchunk, wGuildC, wGuildCompType, wGuildCompPairwise, ...
            wGuildWeightParams, wGuildCompSymmetry, wGuildCompCorr );
        
        
        % CREATE JACOBIANS  -------------------------------------------------------
        
        [ J, repMat ] = ...
            msf_create_Jacobian_ssmnw( ...
            JacobianSum, guildStruct, intra, S, ...
            PPchunk, AAchunk, APchunk, PAchunk, ...
            nlinksBtwG, wGuildC, wGuildCompPairwise, ...
            repMat, t, saveJac);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%   CHECK STABILITY    %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [ resMat ] = ...
            msf_analyze_Jacobian_ssmnw( ...
            J, resMat, t );

    end % of t = replicates for-loop
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SAVING    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    % Saving the Jacobians
    if ( saveJac > 0 )
        save( JacobianFileName, 'repMat')
    end
    
    % Saving resMat, the matrix with the stability results.
    save( resMatFileName, 'resMat')
    
    % Saving workspace parameters - simply good practice in case of debug needs
    clear resMat repMat J
    save( workspaceFileName )   
    
else % Scenario parameter combination is deemed unvalid, resMat will be NaNs.
    
    % Saving resMat, the matrix with the stability results. All NaNs here.
    save( resMatFileName, 'resMat')
    
    clear resMat repMat
    save( workspaceFileName )
    
end % of if valid statement

duration = toc;
fprintf('\n\n>>> TIME REQUIRED FOR SCENARIO NR %i WAS %d<<<\n', row, duration)


end % of function