function [ resMat ] = msf_analyze_Jacobian_ssmnw( J, resMat, t )%( J, resMat, t, repNum )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% if ~isempty(repNum)                                                         % Cluster single job setup.
%     t = 1;
% end

[ stable, ~ ] = obj_calculate_local_stability( J );

resMat(t,1) = stable;

if ( stable)
    [ resilience ] = obj_calculate_resilience( J );
    [ ~, reactivity, ~] = obj_calculate_reactivity( J );
    
    resMat(t,2) = reactivity;
    resMat(t,3) = resilience;
    
else                                                                 % If not locally stable,
    resMat(t,2) = NaN;
    resMat(t,3) = NaN;
end



end

