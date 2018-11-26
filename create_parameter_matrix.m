function [ params ] = create_parameter_matrix(SAVE_DIR, saveName, names, varargin)
%CREATE_PARAMETER_MATRIX creates a fully factorial parameter matrix of the
%variables entered to the function. The parameter matrix will be saved as a
%structure with additional data (variables names and column nr for each
%variable).
%   The parameter matrix will contain the combinations of variable values
%   that will constitute the scenarios in a simulation. For example, assume
%   you have a food web model and wish to run scenarios with 2 different
%   web sizes and 3 different connectances, a total 6 scenarios. Then,
%   as input variables you give first a string array containing the
%   names of your variables, e.g. {'S','C'}, and then two vectors
%   containing the levels of your variables, e.g. [10 20] for the species
%   richness and [0.2 0.5 0.8] for the connectance. Finally, you put in, as
%   a string, the path to where you want your parametermatrix saved.
%
%   @AUTHORS
%   Alva Curtsdotter, PhD student, IFM, Link?ping University, 2014-01-16

% Check input.
if ( length(names) ~= length(varargin) )
        sprintf('\n>>> ERROR! >>>\n\n >>> In create_parameter_matrix:\n >>> Missmatch between the Number of variables and the Number of variable names in input! >>>')
end

nv = length(varargin);                                                      % Nr of variables.

nl = zeros(1,nv);                                                           % Get nr of levels of each variable.
for i = 1:nv;
   nl(i) = length(varargin{i}); 
end

ns = zeros(1,nv);
ns(1) = nl(1);
for i = 2:nv
    ns(i) = ns(i-1)*nl(i);
end

matrix = zeros(ns(end),nv);
matrix(:,1) = repmat(sort(varargin{1}),1,ns(end)/ns(1));
for i = 2:nv;
   matrix(:,i) = repmat( sort ( repmat ( varargin{i},1,ns(i-1) ) ),1,ns(end)/ns(i)  ); 
end

if size(unique(matrix, 'rows'),1) ~= size(matrix,1)
    sprintf('\n>>> ERROR! >>>\n\n >>> In create_parameter_matrix:\n >>> Creating parameter matrix failed. Identical rows (scenarios) were created. >>>')
end

matrix = [ [1:ns(end)]' matrix];

params.names = names;
params.data = matrix;
for i = 1:length(names)
   params = setfield(params, names{i}, i+1) ;
end

if ~isempty(SAVE_DIR)
    saveName = [SAVE_DIR, saveName];
    save(saveName,'params')
end

end