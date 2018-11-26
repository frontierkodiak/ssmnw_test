function [ PPchunk, AAchunk ] = msf_create_weighted_wGuildStructure_ssmnw( ...
    APchunk, wGuildC, wGuildCompType, pairwise, weightParams,...
    wGuildCompSymmetry, wGuildCompCorr )

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

P = size(APchunk,2);                                                        % Number of species in guild 1 (plants).
A = size(APchunk,1);                                                        % Number of species in guild 2 (animals).

if wGuildC == 0
    
    PPchunk = zeros(P);                                                     % Create plant-plant interaction block for zero connectance case.
    AAchunk = zeros(A);                                                     % Create animal-animal interaction block for zero connectance case.
    
else
    
    % Number of links to distribute
    
    nlinksPGuild    = round(wGuildC*P*P);                                   % Number of links within guilds.
    nlinksPGuild    = min(nlinksPGuild, P*(P-1));                           % Not distributing intraspecific links here. So max nlinks allowed P*(P-1).
    PPchunk         = zeros(P,P);                                           % Initiate withinGuild interaction matrix.
    
    nlinksAGuild    = round(wGuildC*A*A);                                   % Number of links within guilds.
    nlinksAGuild    = min(nlinksAGuild, A*(A-1));                           % Not distributing intraspecific links here. So max nlinks allowed A*(A-1).
    AAchunk         = zeros(A,A);                                           % Initiate withinGuild interaction matrix.
    
    if ( pairwise && round(nlinksPGuild/2) ~= nlinksPGuild/2 )              % If odd number of interactions within guild is required but they need to be pairwise...
        rv = rand(1)-1/2;
        sign = rv/abs(rv);
        nlinksPGuild = nlinksPGuild + sign;
        nlinksAGuild = nlinksAGuild - sign;
        
        if ( nlinksPGuild+nlinksAGuild ~= round(wGuildC*P*P)+round(wGuildC*A*A) )
            error('Error in msf_create_weighted_wGuildStructure_ssmnw: Incorrect number of interactions.')
        end
    end
    
    % Create binary structure
    
    if pairwise                                                             % If binary interactions are to be pairwise.
        
        % Plant binary interacion structure
        allowedIdx = find(tril(ones(P), -1));                               % Get indices for lower triangle elements.
        idx = allowedIdx(randperm(length(allowedIdx), nlinksPGuild/2));     % Randomly draw elements for half of the interactions.
        [r,c] = ind2sub(size(PPchunk),idx);                                 % Get the corresponding upper triangle elements for the pairwise interactions.
        PPchunk(idx) = 1;                                                   % Assign the interactions.
        PPchunk(sub2ind(size(PPchunk), c, r)) = 1;                          % Assign the interactions.
        
        % Bug check (check actual pairwise-ity)
        [ur,uc] = find(triu(PPchunk,1)); urc=[ur,uc]; urc=sortrows(urc);
        [lr,lc] = find(tril(PPchunk,-1)); lcr=[lc,lr]; lcr=sortrows(lcr);
        if ~all(lcr(:)==urc(:))
            fprintf('\n>>> ERROR in msf_create_weighted_wGuildStructure_ssmnw! \n>>> Failure to create pairwise within guild interaction structure!')
            return
        end
        
        % Animal binary interacion structure
        allowedIdx = find(tril(ones(A), -1));                               % Get indices for lower triangle elements.
        idx = allowedIdx(randperm(length(allowedIdx), nlinksAGuild/2));     % Randomly draw elements for half of the interactions.
        [r,c] = ind2sub(size(AAchunk),idx);                                 % Get the corresponding upper triangle elements for the pairwise interactions.
        AAchunk(idx) = 1;                                                   % Assign the interactions.
        AAchunk(sub2ind(size(AAchunk), c, r)) = 1;                          % Assign the interactions.
        
        % Bug check (check actual pairwise-ity)
        [ur,uc] = find(triu(AAchunk,1)); urc=[ur,uc]; urc=sortrows(urc);
        [lr,lc] = find(tril(AAchunk,-1)); lcr=[lc,lr]; lcr=sortrows(lcr);
        if ~all(lcr(:)==urc(:))
            fprintf('\n>>> ERROR in msf_create_weighted_wGuildStructure_ssmnw! \n>>> Failure to create pairwise within guild interaction structure!')
            return
        end
        
    else                                                                    % If interactions do not have to be pairwise.
        fprintf('\n>>> ERROR in msf_create_weighted_wGuildStructure_ssmnw! \n>>> Code cannot (should not) accomodate parameter value of "pairwise" to be zero.')
        % If you allow the non-pairwise option again, be aware that it is only
        % possible with guildStructure, but not without. If you try, the code will
        % ground to a halt in msf_create_jacobian. Develop it if you for some
        % reason want this functionality!
        %         allowedIdx      = find(~eye(P));                                        % Finding index for off-diagonal elements.
        %         PPchunk(allowedIdx(randperm(numel(allowedIdx),nlinksPGuild))) = 1;      % Place links randomly.
        %         allowedIdx      = find(~eye(A));                                        % Finding index for off-diagonal elements.
        %         AAchunk(allowedIdx(randperm(numel(allowedIdx),nlinksAGuild))) = 1;      % Place links randomly.
        
    end
    
    % Create rank-order structure
    
    if ( wGuildCompSymmetry == 2 ||  wGuildCompSymmetry == 3 ) % Prespecified numerical correlation or rank-order correlation
        
        % Create rank-order
        trilPPchunk = tril(PPchunk);
        trilPPchunk(trilPPchunk == 1) = randperm(nlinksPGuild/2);
        triuPPchunk = trilPPchunk';
        
        trilAAchunk = tril(AAchunk);
        trilAAchunk(trilAAchunk == 1) = randperm(nlinksAGuild/2);        % Assign order for weights.
        triuAAchunk = trilAAchunk';
        
    end
    
    % Get weights
    
    if wGuildCompType == 1 % diffuse. This case doesn't make much sense anymore.
        
        Pweights = -weightParams(2)*ones(nlinksPGuild,1);
        Aweights = -weightParams(2)*ones(nlinksAGuild,1);
        
    elseif wGuildCompType == 2 % other than diffuse (used to be diffuseRand)
        
        if wGuildCompSymmetry == 2 % Prespecified numerical correlation
            
            [ Pweights, rhoOutP ] = obj_generate_corrPairs_randVar(...
                'randn', 1, 1, nlinksPGuild/2, wGuildCompCorr,  ...
                repmat(weightParams(1),1,2), repmat(weightParams(2),1,2), 0);
            
            [ Aweights, rhoOutA ] = obj_generate_corrPairs_randVar(...
                'randn', 1, 1, nlinksAGuild/2, wGuildCompCorr,  ...
                repmat(weightParams(1),1,2), repmat(weightParams(2),1,2), 0);
            
            Pweights = -Pweights;
            Aweights = -Aweights;
            
        elseif wGuildCompSymmetry == 3 % Rank-order correlation
        
            Pweights = -sort(abs(weightParams(2).*randn(nlinksPGuild/2,2) + weightParams(1)));
            Aweights = -sort(abs(weightParams(2).*randn(nlinksAGuild/2,2) + weightParams(1)));
            
        else % Independent case or perfect symmetry
            
            Pweights = -abs(weightParams(2).*randn(nlinksPGuild,1) + weightParams(1));
            Aweights = -abs(weightParams(2).*randn(nlinksAGuild,1) + weightParams(1));
            
        end
    end
    
    % Distribute weights
    
    if wGuildCompType == 1 % diffuse. This case doesn't make much sense anymore.
        
        PPchunk(PPchunk==1) = Pweights;
        AAchunk(AAchunk==1) = Aweights;
        
    elseif wGuildCompType == 2 % other than diffuse (used to be diffuseRand)
        
        if ( wGuildCompSymmetry == 0 ||  wGuildCompSymmetry == 1) % Independent case or perfect symmetry case
            
            PPchunk(PPchunk==1) = Pweights;
            AAchunk(AAchunk==1) = Aweights;
            
            if wGuildCompSymmetry == 1 % Perfect symmetry
                
                PPchunk = tril(PPchunk) + tril(PPchunk)';
                AAchunk = tril(AAchunk) + tril(AAchunk)';
                
            end
            
        elseif ( wGuildCompSymmetry == 2 ||  wGuildCompSymmetry == 3 ) % Prespecified numerical correlation or rank-order correlation
            
            % Distribute weights for plant interactions
            for j = 1:nlinksPGuild/2
                trilPPchunk(find(trilPPchunk==j)) = Pweights(j,1);
                triuPPchunk(find(triuPPchunk==j)) = Pweights(j,2);
            end
           
            % Bugcheck of correlation
            lintrilpp = abs(trilPPchunk(trilPPchunk~=0));
            lintriupp = triuPPchunk';
            lintriupp = abs(lintriupp(lintriupp~=0));
            if ( wGuildCompSymmetry == 2 )
                if ( round(corr(lintrilpp, lintriupp),4) ~= round(rhoOutP,4) )
                    error('Error in msf_create_weighted_wGuildStructure_ssmnw: Correlation of pairwise plant interactions is incorrect.')
                end
            elseif ( wGuildCompSymmetry == 3 )
                if ( corr(lintrilpp, lintriupp, 'type', 'Spearman') ~= 1 )
                    error('Error in msf_create_weighted_wGuildStructure_ssmnw: Rank correlation of pairwise plant interactions is incorrect.')
                end
            end
            PPchunk = trilPPchunk + triuPPchunk;
            
            % Distribute weights for animal interactions
            for j = 1:nlinksAGuild/2
                trilAAchunk(find(trilAAchunk==j)) = Aweights(j,1);
                triuAAchunk(find(triuAAchunk==j)) = Aweights(j,2);
            end
            
            % Bugcheck of correlation
            lintrilaa = abs(trilAAchunk(trilAAchunk~=0));
            lintriuaa = triuAAchunk';
            lintriuaa = abs(lintriuaa(lintriuaa~=0));
            if ( wGuildCompSymmetry == 2 )
                if ( round(corr(lintrilaa, lintriuaa),4) ~= round(rhoOutA,4) )
                    error('Error in msf_create_weighted_wGuildStructure_ssmnw: Correlation of pairwise animal interactions is incorrect.')
                end
            elseif ( wGuildCompSymmetry == 3 )
                if ( corr(lintrilaa, lintriuaa, 'type', 'Spearman') ~= 1 )
                    error('Error in msf_create_weighted_wGuildStructure_ssmnw: Rank correlation of pairwise animal interactions is incorrect.')
                end
            end
            AAchunk = trilAAchunk + triuAAchunk;            
            
        end %of wGuildCompSymmetry if statement for weight distribution
        
    end %of wGuildCompType if statement for weight distribution
    
end % of wGuildC if statement

end

