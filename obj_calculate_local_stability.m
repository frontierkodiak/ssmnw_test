function [ out1, out2, out3 ] = obj_calculate_local_stability( J )

%OBJ_CALCULATE_LOCAL_STABILITY calculates local stability for a given symmetric 
% Jacobian matrix, J.
%
%   @INPUT
%   J         A Jacobian matrix. Symmetric M-by-N matrix.
%
%   @OUTPUT
%   out1  Local stability. Logical scalar. 1 if system is locally stable; 0 if not.
%   out2  The real part of all eigenvalues. Vector.
%   out3  All eigenvalues, in full. Vector.
%   
%   @AUTHORS
%   Alva Curtsdotter, Post doc @ BrosiLab, Dep of Environmental Sciences,
%   Emory University, Atlanta, Georgia, USA. Code initiated 2017-11-29.
%
%--------------------------------------------------------------------------

eigVec  = eig(J);                                                               % Get eigenvalues. Vector.
realVec = real(eigVec);                                                         % Extract the real part of each eigenvalue. Vector.
stable  = all(realVec < 0 );                                                    % Check if all eigenvalus have a negative real part.

out1 = stable;
out2 = realVec;
out3 = eigVec;

end % of function 