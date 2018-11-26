function [ J, repMat ] = ...
    msf_create_Jacobian_ssmnw( ...
    JacobianSum, guildStruct, intra, S, PPchunk, AAchunk, APchunk, PAchunk, ...
    nlinksBtwG, wGuildC, pairwise, repMat, t, saveJac )

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

P = size(PPchunk, 1);
A = size(AAchunk, 1);

if JacobianSum == 0 % Sum of Jacobian not forced to be constant across scenarios.
    
    if guildStruct == 0                  % Scramble Jacobian if noBlock case
        
        intra = intra*ones(1,S);
        J     = diag(intra);
        
        PPchunk = triu(PPchunk,1)+tril(PPchunk,-1);
        PPvals = PPchunk(PPchunk~=0);                       % Values of within-guild interactions (plants or guild 1)
        
        AAchunk = triu(AAchunk,1)+tril(AAchunk,-1);
        AAvals = AAchunk(AAchunk~=0);                       % Values of within-guild interactions (animals or guild 2)
        
        [AProws, APcols, APvals] = find(APchunk);           % Values of between-guild interactions (effects of plants on animals)
        apstuff = sortrows([AProws, APcols, APvals]);
        
        [PArows, PAcols, PAvals] = find(PAchunk);           % Values of between-guild interactions (effects of animals on plants)
        pastuff = sortrows([PAcols, PArows, PAvals]);
        
        if ( all(pastuff(:,1)==apstuff(:,1)) && all(pastuff(:,2)==apstuff(:,2)) )            
            logiVec = rand(length(PArows),1);
            logiVec(logiVec>0.5)  = 1;
            logiVec(logiVec<=0.5) = 0;
            btwGuildIS = logiVec.*[pastuff(:,3),apstuff(:,3)] + ~logiVec.*[apstuff(:,3),pastuff(:,3)]; % Values of all between-guild interactions, organized so as to keep interaction type and symmetry intact.
        else
            fprintf('\n>>> ERROR in msf_create_Jacobian_ssmnw! \n>>> Failure to create noBlock structure! Simulation ended @btwGuildIS.')
            return
        end
        
        % Redistribute values of within-guild interactions (plants or guild 1) across entire Jacobian
        if ~isempty(PPvals)
            
            if pairwise
                
                tmpJ = triu(J)+tril(ones(size(J)),-1);
                ze =(find(tmpJ==0));
                idx = ze(randperm(length(ze), length(PPvals)/2));
                [r,c] = ind2sub(size(J),idx);
                J(idx) = PPvals(1:length(PPvals)/2);
                J(sub2ind(size(J), c, r)) = PPvals( (1+length(PPvals)/2) : length(PPvals) );
                
                % Bug check
                [ur,uc] = find(triu(J,1)); urc=[ur,uc]; urc=sortrows(urc);
                [lr,lc] = find(tril(J,-1)); lcr=[lc,lr]; lcr=sortrows(lcr);
                if ~all(lcr(:)==urc(:))
                    fprintf('\n>>> ERROR in msf_create_Jacobian_ssmnw! \n>>> Failure to create noBlock structure! Simulation ended @PPvals.')
                    return
                end
                
            else
                fprintf('\n>>> ERROR in msf_create_Jacobian_ssmnw! \n>>> Code cannot accomodate parameter value of "pairwise" to be zero.')
                clear J
                return
            end
            
        end
        
        % Redistribute values of within-guild interactions (animals or guild 2) across entire Jacobian
        if ~isempty(AAvals)
            
            if pairwise
                
                tmpJ = triu(J)+tril(ones(size(J)),-1);
                ze =(find(tmpJ==0));
                idx = ze(randperm(length(ze), length(AAvals)/2));
                [r,c] = ind2sub(size(J),idx);
                J(idx) = AAvals(1:length(AAvals)/2);
                J(sub2ind(size(J), c, r)) = AAvals( (1+length(AAvals)/2) : length(AAvals) );
                
                % Bug check
                [ur,uc] = find(triu(J,1)); urc=[ur,uc]; urc=sortrows(urc);
                [lr,lc] = find(tril(J,-1)); lcr=[lc,lr]; lcr=sortrows(lcr);
                if ~all(lcr(:)==urc(:))
                    fprintf('\n>>> ERROR in msf_create_Jacobian_ssmnw! \n>>> Failure to create noBlock structure! Simulation ended @AAvals.')
                    return
                end
                
            else
                fprintf('\n>>> ERROR in msf_create_Jacobian_ssmnw! \n>>> Code cannot accomodate parameter value of "pairwise" to be zero.')
                clear J
                return
            end
            
        end
        
        % Redistribute values of between-guild interactions across entire Jacobian, while keeping interaction type and symmetry intact.
        tmpJ = triu(J)+tril(ones(size(J)),-1);
        ze =(find(tmpJ==0));
        idx = ze(randperm(length(ze), nlinksBtwG));
        [r,c] = ind2sub(size(J),idx);        
        J(idx) = btwGuildIS(:,1);
        J(sub2ind(size(J), c, r)) = btwGuildIS(:,2);
        
    else                                                    % Create Jacobian with block structure for Oarrays other than 'noBlock'
        
        if wGuildC == 0
            intra = intra*ones(1,S);
            J     = diag(intra);
        else
            J            = [PPchunk, zeros(P,A); zeros(A,P), AAchunk];
            J(1:S+1:end) = intra;
        end
        
        J(P+1:S,1:P) = APchunk;
        J(1:P,P+1:S) = PAchunk;
        
    end
    
else
    
    fprintf('\n\n>>> ERROR in msf_create_Jacobian_ssmnw: Code not written to accomodate JacobianSum being constant (value %i). <<<\n', JacobianSum)
    return
    
end

% Put the Jacobian in the storage structure
if ( ismember(t, 1:saveJac) )
    repMat(:,:,t) = J;
end

end

