function [ out1, out2 ] = obj_calculate_resilience( J )

%OBJ_CALCULATE_RESILIENCE calculates engineering, or asymptotic, resilience for 
% a given symmetric Jacobian matrix, J, that is already known to be locally stable,
% i.e. all eigenvalues have negative real parts. 
%
%   @INPUT
%   J         A Jacobian matrix. Symmetric M-by-N matrix.
%
%   @OUTPUT
%   out1  Asymptotic resilience, i.e. larger value equals faster return rate to equilibrium equals higher stability. Scalar.
%   out2  Characteristic return time. Inverse of the above. Higher value equals longer return time equals lower stability. Scalar.
%   
%   @AUTHORS
%   Alva Curtsdotter, Post doc @ BrosiLab, Dep of Environmental Sciences,
%   Emory University, Atlanta, Georgia, USA. Code initiated 2017-11-30.
%
%--------------------------------------------------------------------------

eigVec   = eig(J);                                                                % Get eigenvalues. Vector.
realVec  = real(eigVec);                                                          % Extract the real part of each eigenvalue. Vector.
largest  = max(realVec);                                                          % Get maximum value of real parts. Scalar. 

if ( ~all(realVec < 0 ) )                                                         % Dummy check that Jacobian is locally stable.
  fprintf('\n>>> NOTICE! >>> In obj_calculate_resilience:\n\n\t>>> Provided Jacobian is not locally stable. >>>\n\n' )
  return
end 
%    
% if  ( sum( abs(realVec)==abs(largest) ) > 1 )                                    % Notify user when multiple dominant eigenvalus exist. 
%   fprintf('\n>>> NOTICE! >>> In obj_calculate_resilience:\n\n\t>>> %i dominant eigenvalues. What does this mean for resilience? I currently assume that it matters not for resilience whether there is one or more dominant eigenvalues. >>>\n\n', sum( abs(realVec)==abs(largest) ) )
% end 

resilience = abs(largest);                                                       % Resilience is the absolute value of the largest (i.e. least negative) real part of the eigenvalues. Based on communication w P. STaniczenko, and Arnoldi et al 2016 JTB.

out1 = resilience;                                                                % Assign out1.
out2 = 1/resilience;                                                              % Calculate return time and assign out2.

end % of function 