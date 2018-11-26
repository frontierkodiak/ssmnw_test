function [ APchunk, PAchunk ] = msf_create_weighted_btwGuildStructure_ssmnw( ...
    APchunk, PAchunk, weightParams, nlinksBtwG, numIntTypes, intTypeSymmetry, ...
    intTypeCorr, intTypes )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if numIntTypes == 1         % This is the only case defined so far.
    
    if ( intTypeSymmetry == 0 || intTypeSymmetry == 3 )% Independence of pairwise interaction strengths or Rank-Order Correlation.
        
        APweights = sort(abs(weightParams(2).*randn(nlinksBtwG,1) + weightParams(1))); % Independently drawn weights, sorted by absolute value.
        PAweights = sort(abs(weightParams(2).*randn(nlinksBtwG,1) + weightParams(1)));
                
    elseif intTypeSymmetry == 1  % Perfect symmetry of pairwise interaction strengths
        
        APweights = sort(abs(weightParams(2).*randn(nlinksBtwG,1) + weightParams(1)));
        PAweights = APweights;
        
    elseif intTypeSymmetry == 2  % Prespecified Pearsons correlation of pairwise IS.
        
        [ weights, ~ ] = obj_generate_corrPairs_randVar(...
            'randn', 1, 1, nlinksBtwG, intTypeCorr,  ...
            repmat(weightParams(1),1,2), repmat(weightParams(2),1,2), 0);
        
        APweights = weights(:,1);
        PAweights = weights(:,2);
        
        % Note that this option only makes sense for a random
        % interaction structure (even if you sort, cause you can only sort one
        % column. I.e. you can either pre-specify the rank-order structure of
        % BOTH blocks, OR the pairwise correlations, but not both as far as I
        % can see atm!).        
        
    end % of intTypeSymmetry if statement
    
    if intTypes == 1 %competition
        
        for j = 1:nlinksBtwG
            APchunk(find(APchunk==j)) = -APweights(j);     % Distribute weights according to order.
            PAchunk(find(PAchunk==j)) = -PAweights(j);
        end
        
    elseif intTypes == 2 % mutualism
        
        for j = 1:nlinksBtwG
            APchunk(find(APchunk==j)) = APweights(j);     % Distribute weights according to order.
            PAchunk(find(PAchunk==j)) = PAweights(j);
        end
        
    elseif intTypes == 3 %trophic
        
        for j = 1:nlinksBtwG
            APchunk(find(APchunk==j)) = APweights(j);     % Distribute weights according to order.
            PAchunk(find(PAchunk==j)) = -PAweights(j);
        end
        
    end % of intTypes if statement
    
else
    
    fprintf('\n\n>>> ERROR in msf_create_weighted_btwGuildStructure_ssmnw: Code not written to accomodate %i interaction types! <<<\n', numIntTypes)
    return
    
end % of numIntTypes if statement

end

