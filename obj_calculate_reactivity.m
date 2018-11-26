function [ out1, out2, out3, out4, out5 ] = obj_calculate_reactivity( J )

%OBJ_CALCULATE_REACTIVITY calculates reactivity and initial resilience, i.e. the opposite of reactivity, 
% for a given symmetric Jacobian matrix, J, that is already known to be locally stable,
% i.e. all eigenvalues have negative real parts. Based on Tang & Allesina 2014 Frontiers and 
% Arnoldi et al. 2016 Journal of Theoretical Biology.
%
%   @INPUT
%   J         A Jacobian matrix. Symmetric M-by-N matrix.
%
%   @OUTPUT
%   out1  Whether system is reactive (1) or not (0). Scalar.
%   out2  Reactivity. Scalar.
%   out3  Initial resilience. Scalar.
%   out4  The real part of all eigenvalues. Vector.
%   out5  All eigenvalues, in full. Vector.  
%   
%   @AUTHORS
%   Alva Curtsdotter, Post doc @ BrosiLab, Dep of Environmental Sciences,
%   Emory University, Atlanta, Georgia, USA. Code initiated 2017-11-30.
%
%--------------------------------------------------------------------------

eigVec   = eig(J);                                                                % Get eigenvalues. Vector.
realVec  = real(eigVec);                                                          % Extract the real part of each eigenvalue. Vector.
   
if ( ~all(realVec < 0 ) )                                                         % Dummy check that Jacobian is locally stable.
  fprintf('\n>>> NOTICE! >>> In obj_calculate_reactivity:\n\n\t>>> Provided Jacobian is not locally stable. >>>\n\n' )
  return
end 

eigVec   = eig((J+J')/2);                                                         % Get eigenvalues. Vector.
realVec  = real(eigVec);                                                          % Extract the real part of each eigenvalue. Vector.

reactivity =  max(realVec);                                                       % Calculate reactivity. Scalar. Equation from Tang & Allesina 2014 Frontiers. 
reactive   = reactivity > 0;                                                      % Whether system is reactive or not.

out1 = reactive;
out2 = reactivity;                                                           
out3 = -reactivity;                                                               % Initial resilience. Arnoldi et al. 2016 Journal of Theoretical Biology.
out4 = realVec;                                                               
out5 = eigVec;

end % of function 