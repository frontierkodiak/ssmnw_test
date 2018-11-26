function [ pairs, cf ] = obj_generate_corrPairs_randVar(...
    dist, doAbsolute, doSort, n, rho, my, sigma, doFigs)

%OBJ_GENERATE_CORRPAIRS_ABSVALS_RANDVAR generates two sets of random values
%from a specified distribution, such that the values come in pairs with a
%specified correlation. The paired values must come from the same type of
%of distribution, but can come from distributions with different mean and
%standard deviation.

%   @INPUT
%   dist        The type of distribution to use. String. Accepted values are
%               currently 'randn' and 'rand'.
%   doAbsolute  Take the absolute value of the random variable. Logical: 1=yes; 0=no.
%   doSort      Return pairs sorted by strength of one of the two sets of values 
%               (i.e. one "side" of the pairs). Logical: 1=yes; 0=no.
%   n           The number of pairs. Integer.
%   rho         The desired correlation of paired values. Integer.
%   my          The means of the distributions. Vector [my1 my2].
%   sigma       The standard deviations of the distributions. Vector [sigma1 sigma2]
%   doFigs      Whether to make figure of random variabled (development purposes). Logical: 1=yes; 0=no.
%   
%   @OUTPUT
%   pairs   The generated values. n-by-2 matrix, i.e. n pairs of values.
%
%   @AUTHORS
%   Alva Curtsdotter, Post doc, Department of Environmental Sciences, Emory
%   University, Atlanta, 2018-10-24.
%
%--------------------------------------------------------------------------

% Definitions -------------------------------------------------------------

rhoMatchSensitivity = 2; % Number of digits to which to match target and actual rho.
turnMax  = 10000;        % How many times to try to get reach above criteria before giving up.

% Correlation -------------------------------------------------------------

rhoMatch = 0;
turn     = 1;

while rhoMatch == 0
    
    % Draw 2 sets of random values
    if doAbsolute 
        X = abs(eval( [dist, '(n,2)'] ));
    else
        X = eval( [dist, '(n,2)'] );
    end    
    if doSort
        X = sortrows(X);
    end
    
    % Create a third set of values with pairwise correlations, rho, with set 1.
    X(:,3) = rho.*X(:,1)+sqrt(1-rho^2).*X(:,2);
    
    % Attempt to ensure sufficient match between target and actual rho.
    cf  = corr(X(:,1),X(:,3));
    cfr = round(cf, rhoMatchSensitivity);    
    if ( cfr == rho || turn == turnMax )
        rhoMatch = 1;
        if ( turn == turnMax )
            fprintf('\n\n>>> ERROR: In obj_correlated_pairs...: Requested (%d) and actual (%d) correlation do not match! <<<\n', rho, cfr)
            fprintf('\n\n>>> WARNING: In obj_correlated_pairs...: CONTINUES anyway. Reached maximum number of trials: %i. <<<\n', turnMax )
        end
    else
        fprintf('\n\n>>> In obj_correlated_pairs...: Requested (%d) and actual (%d) correlation do not match! Turn %i of while-loop. <<<\n', rho, cfr, turn)
        turn     = turn+1;
    end    
end


% Transformation ----------------------------------------------------------

% Transform from standard distribution to the specified means and standard
% deviations. This does not effect the pairwise correlation.
Y1 = (my(1) + sigma(1).*X(:,1));
Y2 = (my(2) + sigma(2).*X(:,3));

pairs = [Y1, Y2];


% Analysis (developmental) ------------------------------------------------

% If you want to do visual investigation of values
if doFigs
    figure()
    plot(X(:,1),X(:,2), 'g.')
    
    figure()
    plot(X(:,1),X(:,3), 'g.')
    
    figure()
    plot(Y1,Y2, 'b.')
    cfy = corr(Y1,Y2);
    disp(cfy)
end


end

