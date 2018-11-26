function [ APchunk, PAchunk, nlinksBtwG ] = ...
    msf_create_binary_btwGuildStructure_ssmnw( project, S, numGuilds, btwGuildC, ...
    intTypeSymmetry, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if numGuilds == 2 % This is the only case defined so far.
    
    % Define basic network parameters
    A     = round(S/2);                   % Number of species in Guild 1.
    P     = S - A;                        % Number of species in Guild 2.
    APchunk    = zeros(A,P);              % Initiate incidence matrix.
    nlinksBtwG = round(btwGuildC*A*P);    % Number of links between guilds.
    
    % Doublecheck that matrix is symmetric (an assumption of the function)
    if (A ~= P )
        fprintf('\n>>> ERROR! >>>\n\n >>> In msf_create_binary_btwGuildStructure_ssmnw:\n >>> Lacking functionality for different sized guilds. >>>\n')
        return
    end
    
    % Project 1: guild structure  
    if isempty(varargin)
        
        % Dummy check, that we want to run project 1.
        if project ~= 1
            fprintf('\n\n>>> ERROR in msf_create_binary_btwGuildStructure_ssmnw: Wrong number of input variables for project %i! <<<\n', project)
            return
        end
        
        % Create binary random interaction structure.
        APchunk(sub2ind(size(APchunk), randperm(A), randperm(P))) = 1;           % Make sure each species has a btwGuild interaction.
        idx = find(APchunk == 0);                                                % Find indices of the remaining elements (i.e. which elements that have not already been assigned a 1 for interation)
        APchunk(idx(randperm(length(idx),nlinksBtwG-sum(sum(APchunk))))) = 1;    % Assign the remaining number of interactions.
        PAchunk = APchunk';
        
        % Create rank-order
        APchunk(APchunk == 1) = randperm(nlinksBtwG);                            % Assign order for weights.
        if intTypeSymmetry == 0
            PAchunk(PAchunk == 1) = randperm(nlinksBtwG);                        % Assign order for weights; order in one block independent of that in another.
        else
            PAchunk = APchunk';                                                  % Rank-order correlation = 1. necessary/convenient also for creating prespecified numerical correlations other than zero.
        end
        
    else % Project 2: nestedness 
        
        % Dummy check, that we want to run project 2.
        if project ~= 2
            fprintf('\n\n>>> ERROR in msf_create_binary_btwGuildStructure_ssmnw: Wrong number of input variables for project %i! <<<\n', project)
            return
        end
        
        % Define additional basic network parameters, specific to project 2.
        bNest = varargin{1};        % Level of binary nestedness.
        qNest = varargin{2};        % Level of quantitative nestedness.
        
        if bNest == 1 % Perfect binary nestedness.
            
            APchunk = zeros(A, P); 
            APchunk(:, P) = 1;
            APchunk(1,:) = 1;
            nlinksRemain = nlinksBtwG-sum(sum(APchunk));
            
            subAPchunk = zeros(A-1, P-1); 
            kVec = [(A-2):-1:0 -(1:(A-2))];       % Vector of diagonal number (we want to loop over all the diagonals of the fillOrderMatrix)
            
            for kidx = 1:length(kVec)                       % Loop over all the diagonals.
                k = kVec(kidx);                             % k is the diagonal we're currently working with.
                diagIdx = false(A-1-abs(k),1);            % Vector of ones, of the same length as diagonal k.
                diagIdx(1:min(length(diagIdx), nlinksRemain)) = 1;            % Vector of ones, of the same length as diagonal k.
                               
                subAPchunk(diag(diagIdx,k)) = 1;       % Set the elements of the current diagonal...
                    % ...to the correct elements in orderElemtens.
                nlinksRemain = nlinksRemain-length(diagIdx);
                if nlinksRemain <= 0
                    break 
                end
            end
            
            APchunk(~APchunk) = subAPchunk;
            APchunk = fliplr(APchunk);         
                       
        elseif bNest == 0 % Random binary network.
            
            APchunk(sub2ind(size(APchunk), randperm(A), randperm(P))) = 1;
            idx = find(APchunk == 0);
            APchunk(idx(randperm(length(idx),nlinksBtwG-sum(sum(APchunk==1))))) = 1;         % Place links randomly.
            
        elseif bNest == -1 % Perfectly antinested network.
            
%             A = 8; nSinG2 = A;
%              APchunk    = zeros(A,nSinG2);              % Initiate incidence matrix.
%             nlinksBtwG = 32;    % Number of links between guilds.
%     
            meanDegree = nlinksBtwG/(P);
            degreeVec = floor(meanDegree)*ones(1,P);
            diff = nlinksBtwG-sum(degreeVec);
            degreeVec(1:diff) = degreeVec(1:diff)+1;
            degreeVec(P+1) = 0; %dummy
            
            rf = 1; % first row to fill for c  
            for c = 1:P
                ridx = rf:(rf+degreeVec(c)-1);
                foo  = ridx-P;
                ridx(foo > 0) = foo(foo > 0);
                APchunk(ridx,c) = 1;
                rf = (rf+degreeVec(c)-1) - (degreeVec(c+1)-1) +1; %for next c                
            end
            APchunk = fliplr(APchunk);
%             APchunk
%             sum(APchunk,1)
%             sum(APchunk,2)
            
        else
            fprintf('\n\n>>> ERROR in msf_create_binary_btwGuildStructure_ssmnw: Lacking functionality for this level of binary nestedness: %i! <<<\n', bNest)
            return
        end
        
        % Sort it so as not to underestimate its nestedness (important for
        % when adding the weights). Should be redundant for bNest = 1 but
        % necessary otherwise.
        [~,colIdx] = sort(sum(APchunk,1), 'descend');
        [~,rowIdx] = sort(sum(APchunk,2), 'descend');
        APchunk = APchunk(:, colIdx);
        APchunk = APchunk(rowIdx, :);        
        
        if sum(sum(APchunk))~=nlinksBtwG
               fprintf('\n\n>>> ERROR in msf_create_binary_btwGuildStructure_ssmnw: Wrong number of links in APchunk! <<<\n')
        end
            
        if qNest == 1 % Perfect quantitative nestedness.
            
            fillOrderMatrix = zeros(A, P);
            orderElements = numel(fillOrderMatrix):-1:1;
            oECounter = 1;
            kVec = [(A-1):-1:0 -(1:(A-1))];
            for kidx = 1:length(kVec)
                k = kVec(kidx);
                diagIdx = true(A-abs(k),1);
                fillOrderMatrix(diag(diagIdx,k)) = orderElements(oECounter:oECounter+length(diagIdx)-1);
                oECounter = oECounter+length(diagIdx);
            end
            fillOrderMatrix = fliplr(fillOrderMatrix);
            
        elseif qNest == 0 % Random quantitative structure.
            
            fillOrderMatrix = reshape(randperm(A*P),[A, P]);
            
        elseif qNest == -1 % Perfect quantitative antinestedness.
            
            fillOrderMatrix = zeros(A, P);
            orderElements = 1:1:numel(fillOrderMatrix);
            oECounter = 1;
            kVec = [(A-1):-1:0 -(1:(A-1))];
            for kidx = 1:length(kVec)
                k = kVec(kidx);
                diagIdx = true(A-abs(k),1);
                fillOrderMatrix(diag(diagIdx,k)) = orderElements(oECounter:oECounter+length(diagIdx)-1);
                oECounter = oECounter+length(diagIdx);
            end
            fillOrderMatrix = fliplr(fillOrderMatrix);
            
        else
            fprintf('\n\n>>> ERROR in msf_create_binary_btwGuildStructure_ssmnw: Lacking functionality for this level of quantitative nestedness: %i! <<<\n', qNest)
            return
        end
        
        APchunk(APchunk==1) = fillOrderMatrix(APchunk==1);
        [~,idx] = maxk(APchunk(:),nlinksBtwG);
        APchunk(idx) = nlinksBtwG:-1:1;
        
        PAchunk = APchunk';
        
    end % of isempty(varargin)
    
else
    fprintf('\n\n>>> ERROR in msf_create_binary_btwGuildStructure_ssmnw: Code not written to accomodate %i guilds<<<\n', numGuilds)
    return
end % of if numGuilds==2

end %of function

